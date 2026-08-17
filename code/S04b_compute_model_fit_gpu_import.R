# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# Define all strategies required for diagnosis (matches S04a)
validation_strategies <- c("metso_holdout", "north_south") # "route_blocked_cv") #, "random_cv")

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")
run_parallel <- TRUE

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 4. Main Loop: Iterate Over MCMC config and Strategies ---
message(paste0("--- Import for Model ID: ", run_config$model_id))

for(i in 1:nrow(mcmc_params)){
  # i <- 1
  run_name <- generate_run_name(run_config)[i]
  base_model_name <- run_name 
  mcmc_params_i <- mcmc_params[i, ]
  
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  fitted_model_path <- file.path("models", base_model_name, paste0("fitted_", base_model_name, ".rds"))
  if(!file.exists(fitted_model_path)) {
    message("  Full fitted model not found. Skipping to next MCMC config.")
    next
  }
  hM_full <- readRDS(fitted_model_path)
  
  hM_full <- Hmsc::alignPosterior(hM_full)
  
  dir_fit <- file.path(models_dir, base_model_name, "model_fit") 
  dir.create(dir_fit, showWarnings = FALSE, recursive = TRUE)
  
  # --- Explanatory Power and WAIC Checks ---
  waic_output_path <- file.path(dir_fit, paste0("waic_", base_model_name, ".rds"))
  expl_output_path <- file.path(dir_fit, paste0("mf_", base_model_name, ".rds"))
  
  if (!file.exists(waic_output_path) || !file.exists(expl_output_path)) {
    message("  Calculating WAIC and Explanatory Model Fit (Full Model)...")
    
    if (!file.exists(waic_output_path)) {
      WAIC_full <- computeWAIC(hM_full)
      saveRDS(WAIC_full, file = waic_output_path)
      message("    WAIC saved to: ", waic_output_path)
      rm(WAIC_full)
    } else {
      message("    WAIC file already exists. Skipping calculation.")
    }
    
    if (!file.exists(expl_output_path)) {
      predY_full <- predict(hM_full, expected = TRUE)
      predY_full <- simplify2array(predY_full)
      MF_explanatory <- Hmsc::evaluateModelFit(hM = hM_full, predY = predY_full)
      saveRDS(MF_explanatory, file = expl_output_path)
      message("    Explanatory Model Fit saved to: ", expl_output_path)
      rm(predY_full, MF_explanatory)
    } else {
      message("    Explanatory Model Fit file already exists. Skipping calculation.")
    }
    gc()
  } else {
    message("  WAIC and Explanatory Model Fit files already exist. Skipping full model processing.")
  }
  
  #--- cross validation evaluation
  for(strategy in validation_strategies) {
    # strategy <- "route_blocked_cv"
    message(paste0("\n  Processing strategy: ", strategy))
    
    # 4.1 Assign label
    if(strategy == "metso_holdout"){
      label <- "ho_metso"
    } else if (strategy == "route_blocked_cv") {
      label <- "cv_route"
    } else if (strategy == "random_cv") {
      label <- "cv_random"
    } else if (strategy == "north_south") {
      label <- "north_south"
    }
    
    # Check if cross-validation evaluation output file already exists
    mf_output_path <- file.path(dir_fit, paste0("mfeval_", base_model_name, "_", label, "_", ".rds"))
    if (file.exists(mf_output_path)) {
      message("    Evaluation file already exists for ", label, ". Skipping strategy execution loop.")
      next
    }
    
    eval_dir <- file.path(models_dir, base_model_name, label)
    
    # 4.2 Load the partition generated in S04a
    partition_path <- file.path(eval_dir, "partition.rds")
    if(!file.exists(partition_path)){
      message("  Partition file not found for ", label, ". Ensure S04a ran successfully. Skipping.")
      next
    }
    partition <- readRDS(partition_path)
    parts <- sort(unique(partition))
    
    # Initialize the prediction array dynamically
    predY <- NULL
    import_success <- TRUE
    
    # 4.3 Combined Import and Predict Loop
    for(p in 1:length(parts)){
      # p <- 2
      if(strategy == "metso_holdout" | strategy == "north_south" ){
        if(p != 1){
          next
        }
      }
      
      message("    Importing and predicting fold: ", parts[p])
      
      res_path <- import_posterior(
        mcmc = run_config$cv$mcmc, 
        config = run_config, 
        run_nm = base_model_name, 
        partition_number = parts[p], 
        label = label
      )
      
      if(is.logical(res_path) && res_path == FALSE) {
        warning("    Failed to import fold ", parts[p], ". Check Python logs. Skipping strategy.")
        import_success <- FALSE
        break
      } 
      
      m_fold <- readRDS(res_path)
      m_fold <- Hmsc::alignPosterior(m_fold)
      postList_fold <- Hmsc::poolMcmcChains(m_fold$postList, start = 1, thin = 1)
      
      # Initialize predY on the first successful fold import
      if(is.null(predY)) {
        postN <- length(postList_fold)
        predY <- array(NA, dim = c(hM_full$ny, hM_full$ns, postN))
      }
      
      if(strategy == "metso_holdout" | strategy == "north_south" ){
        val_idx <- partition == 2 
      }else{
        val_idx <- partition == parts[p] 
      }
      
      
      XData_val <- hM_full$XData[val_idx, , drop = FALSE]
      dfPi_val <- droplevels(hM_full$dfPi[val_idx, , drop = FALSE])
      
      # Execution routing block for predictive processing
      if (run_parallel && parallel::detectCores() > 1) {
        allocated_cores <- max(1, floor(parallel::detectCores() / 2))
        num_cores <- min(allocated_cores, length(postList_fold))
        
        post_chunks <- split(postList_fold, cut(seq_along(postList_fold), num_cores, labels = FALSE))
        
        message("    Predicting across ", num_cores, " cores...")
        
        # Initialize background cluster sockets to prevent GPU/reticulate fork corruption
        cl <- parallel::makeCluster(num_cores)
        
        # Export 
        parallel::clusterExport(cl, 
                                varlist = c("m_fold", "XData_val", "dfPi_val"), 
                                envir = environment())
        
        parallel::clusterEvalQ(cl, library(Hmsc))
        
        pred_fold_chunks <- parallel::parLapply(cl, post_chunks, function(single_chunk) {
          predict(object = m_fold, 
                  post = single_chunk,
                  XData = XData_val, 
                  studyDesign = dfPi_val, 
                  mcmcStep = 1, 
                  expected = TRUE,
                  nParallel = 1)
        })
        
        parallel::stopCluster(cl)
        
        if (any(sapply(pred_fold_chunks, inherits, "try-error"))) {
          stop("Parallel prediction failed in one or more cluster workers.")
        }
        
        pred_fold_combined <- do.call(c, pred_fold_chunks)
        predY[val_idx, , ] <- simplify2array(pred_fold_combined)
        
        rm(pred_fold_chunks, pred_fold_combined, post_chunks)
        
      } else {
        message("    Predicting sequentially on a single core processing pipeline...")
        
        pred_fold <- predict(m_fold, 
                             post = postList_fold, 
                             XData = XData_val, 
                             studyDesign = dfPi_val, 
                             mcmcStep = 1, 
                             expected = TRUE)
        
        predY[val_idx, ,] <- simplify2array(pred_fold)
        
        rm(pred_fold)
      }
      
      # memory cleanup per fold
      rm(m_fold, postList_fold, XData_val, dfPi_val)
      gc()
    }
    
    if(!import_success) next
    
    # --- 4.4 Evaluate Model Fit ---
    message("    Evaluating Model Fit for ", label)
    MFEVAL <- evaluateModelFitRobust(hM = hM_full, predY = predY)
    
    # Save the resulting Evaluation Object
    saveRDS(MFEVAL, file = mf_output_path)
    message("    Successfully saved evaluation to: ", mf_output_path)
    
    # Clean memory before transitioning to the next strategy
    rm(predY, MFEVAL)
    if(exists("predY_auc")) rm(predY_auc)
    gc()
  }
}

message("    End time: ", timestamp())