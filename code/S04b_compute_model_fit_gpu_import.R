# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# Define all strategies required for diagnosis (matches S04a)
validation_strategies <- c( "route_blocked_cv", "random_cv", "metso_holdout") #c("metso_holdout")# 

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")

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
  
  for(strategy in validation_strategies) {
    # strategy <- c("metso_holdout")
    message(paste0("\n  Processing strategy: ", strategy))
    
    # 4.1 Assign label
    if(strategy == "metso_holdout"){
      label <- "ho_metso"
    } else if (strategy == "route_blocked_cv") {
      label <- "cv_route"
    } else if (strategy == "random_cv") {
      label <- "cv_random"
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
      # p <- 1
      
      if(strategy == "metso_holdout" && p != 1) {
        break
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
      
      # Define validation indices based on strategy
      #if (strategy == "metso_holdout") {
      #  val_idx <- partition == 2 
      #} else {
        val_idx <- partition == parts[p] 
      #}
      
      XData_val <- hM_full$XData[val_idx, , drop = FALSE]
      dfPi_val <- droplevels(hM_full$dfPi[val_idx, , drop = FALSE])
      
      pred_fold <- predict(m_fold, 
                           post = postList_fold, 
                           XData = XData_val, 
                           studyDesign = dfPi_val, 
                           mcmcStep = 1, 
                           expected = TRUE)
      
      predY[val_idx, , ] <- simplify2array(pred_fold)
      
      # memory cleanup per fold
      rm(m_fold, postList_fold, pred_fold, XData_val, dfPi_val)
      gc()
    }
    
    if(!import_success) next
    
    # 4.4 Evaluate Model Fit
    message("    Evaluating Model Fit for ", label)
    MFEVAL <- Hmsc::evaluateModelFit(hM = hM_full, predY = predY)
    
    mf_output_path <- file.path(eval_dir, paste0("MF_", base_model_name, "_", label, "_",".rds"))
    saveRDS(MFEVAL, file = mf_output_path)
    message("    Successfully saved evaluation to: ", mf_output_path)
    
    rm(predY, MFEVAL)
    gc()
  }
}
