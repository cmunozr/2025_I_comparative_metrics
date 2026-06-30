# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# Define all strategies required for diagnosis
validation_strategies <- c("route_blocked_cv", "metso_holdout", "random_cv")

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")
run_parallel <- TRUE

# Define models
model_id_pa <- "fbs_M015PA"
model_id_abund <- "fbs_M015"

# Reconstruct storage convention
run_name_pa <- paste0(model_id_pa, "_thin_100_samples_1000_chains_4")
run_name_abund <- paste0(model_id_abund, "_thin_150_samples_1000_chains_4")

# Load full abundance model to be used the evaluation reference
fitted_abund_path <- file.path("models", run_name_abund, paste0("fitted_", run_name_abund, ".rds"))
if(!file.exists(fitted_abund_path)) {
  stop("Full fitted abundance model not found. Cannot proceed with evaluation reference.")
}
hM_full_abund <- readRDS(fitted_abund_path)
hM_full_abund <- Hmsc::alignPosterior(hM_full_abund)

# MCMC parameter specs for the import function
mcmc_spec_pa <- data.frame(samples = 1000, thin = 100, n_chains = 4, transient_proportion = 0.5, adapt_nf_proportion = 0.4)
mcmc_spec_abund <- data.frame(samples = 1000, thin = 150, n_chains = 4, transient_proportion = 0.5, adapt_nf_proportion = 0.4)

dir_fit_hurdle <- file.path(models_dir, "fbs_M014_hurdle", "model_fit")
dir.create(dir_fit_hurdle, showWarnings = FALSE, recursive = TRUE)

# --- 3. Main Loop: Iterate Over Strategies ---
for(strategy in validation_strategies) {
  message(paste0("\nProcessing strategy: ", strategy))
  # strategy <- "metso_holdout"
  if(strategy == "metso_holdout"){
    label <- "ho_metso"
  } else if (strategy == "route_blocked_cv") {
    label <- "cv_route"
  } else if (strategy == "random_cv") {
    label <- "cv_random"
  }
  
  mf_output_path <- file.path(dir_fit_hurdle, paste0("mfeval_hurdle_", label, ".rds"))
  if (file.exists(mf_output_path)) {
    message("    Hurdle evaluation file already exists for ", label, ". Skipping strategy.")
    next
  }
  
  # Use the partition vector from the abundance model directory
  eval_dir <- file.path(models_dir, run_name_abund, label)
  partition_path <- file.path(eval_dir, "partition.rds")
  if(!file.exists(partition_path)){
    message("    Partition file not found for ", label, ". Skipping.")
    next
  }
  partition <- readRDS(partition_path)
  parts <- sort(unique(partition))
  
  predY <- NULL
  import_success <- TRUE
  
  # --- 4. Combined Import and Predict Loop per Fold ---
  for(p in 1:length(parts)){
    # p <- 1
    if(strategy == "metso_holdout" && p != 1) next
    
    message("    Importing and predicting fold: ", parts[p])
    
    # Import Presence-Absence Fold Posterior
    res_path_pa <- import_posterior(
      mcmc = mcmc_spec_pa,
      config = run_config,
      run_nm = run_name_pa,
      partition_number = parts[p],
      label = label
    )
    
    # Import Abundance Fold Posterior
    res_path_abund <- import_posterior(
      mcmc = mcmc_spec_abund,
      config = run_config,
      run_nm = run_name_abund,
      partition_number = parts[p],
      label = label
    )
    
    if((is.logical(res_path_pa) && res_path_pa == FALSE) || (is.logical(res_path_abund) && res_path_abund == FALSE)) {
      warning("    Failed to import posteriors for fold ", parts[p], ". Skipping strategy.")
      import_success <- FALSE
      break
    }
    
    m_fold_pa <- readRDS(res_path_pa)
    m_fold_pa <- Hmsc::alignPosterior(m_fold_pa)
    postList_fold_pa <- Hmsc::poolMcmcChains(m_fold_pa$postList, start = 1, thin = 1)
    
    m_fold_abund <- readRDS(res_path_abund)
    m_fold_abund <- Hmsc::alignPosterior(m_fold_abund)
    postList_fold_abund <- Hmsc::poolMcmcChains(m_fold_abund$postList, start = 1, thin = 1)
    
    if(length(postList_fold_pa) != length(postList_fold_abund)) {
      stop("Mismatched posterior sample dimensions between PA and Abundance models.")
    }
    
    # Initialize the prediction array on the first successful fold import
    if(is.null(predY)) {
      postN <- length(postList_fold_pa)
      predY <- array(NA, dim = c(hM_full_abund$ny, hM_full_abund$ns, postN))
    }
    
    val_idx <- partition == parts[p]
    XData_val <- hM_full_abund$XData[val_idx, , drop = FALSE]
    dfPi_val <- droplevels(hM_full_abund$dfPi[val_idx, , drop = FALSE])
    
    # --- 5. Execution Routing: Parallel vs Sequential ---
    if (run_parallel && parallel::detectCores() > 1) {
      allocated_cores <- max(1, floor(parallel::detectCores() / 2))
      num_cores <- min(allocated_cores, length(postList_fold_pa))
      
      post_chunks_pa <- split(postList_fold_pa, cut(seq_along(postList_fold_pa), num_cores, labels = FALSE))
      post_chunks_abund <- split(postList_fold_abund, cut(seq_along(postList_fold_abund), num_cores, labels = FALSE))
      
      message("    Predicting Hurdle components across ", num_cores, " cores...")
      cl <- parallel::makeCluster(num_cores)
      
      parallel::clusterExport(cl, 
                              varlist = c("m_fold_pa", "m_fold_abund", "XData_val", "dfPi_val", "post_chunks_pa", "post_chunks_abund"), 
                              envir = environment())
      parallel::clusterEvalQ(cl, library(Hmsc))
      
      pred_fold_chunks <- parallel::parLapply(cl, 1:num_cores, function(ch_idx) {
        # Generate expected probabilities from PA component
        p_pa <- predict(object = m_fold_pa, 
                        post = post_chunks_pa[[ch_idx]],
                        XData = XData_val, 
                        studyDesign = dfPi_val, 
                        mcmcStep = 1, 
                        expected = TRUE,
                        nParallel = 1)
        
        # Generate expected values from Abundance component
        p_ab <- predict(object = m_fold_abund, 
                        post = post_chunks_abund[[ch_idx]],
                        XData = XData_val, 
                        studyDesign = dfPi_val, 
                        mcmcStep = 1, 
                        expected = TRUE,
                        nParallel = 1)
        
        # Execute sample-by-sample Hurdle multiplication
        p_hurdle <- lapply(seq_along(p_pa), function(s) {
          p_pa[[s]] * p_ab[[s]]
        })
        return(p_hurdle)
      })
      
      parallel::stopCluster(cl)
      
      if (any(sapply(pred_fold_chunks, inherits, "try-error"))) {
        stop("Parallel prediction failed in one or more cluster workers.")
      }
      
      pred_fold_combined <- do.call(c, pred_fold_chunks)
      predY[val_idx, , ] <- simplify2array(pred_fold_combined)
      rm(pred_fold_chunks, pred_fold_combined, post_chunks_pa, post_chunks_abund)
      
    } else {
      message("    Predicting sequentially on a single core processing pipeline...")
      
      pred_fold_pa <- predict(m_fold_pa, post = postList_fold_pa, XData = XData_val, studyDesign = dfPi_val, mcmcStep = 1, expected = TRUE)
      pred_fold_abund <- predict(m_fold_abund, post = postList_fold_abund, XData = XData_val, studyDesign = dfPi_val, mcmcStep = 1, expected = TRUE)
      
      pred_fold_hurdle <- lapply(seq_along(pred_fold_pa), function(s) {
        pred_fold_pa[[s]] * pred_fold_abund[[s]]
      })
      
      predY[val_idx, , ] <- simplify2array(pred_fold_hurdle)
      rm(pred_fold_pa, pred_fold_abund, pred_fold_hurdle)
    }
    
    rm(m_fold_pa, m_fold_abund, postList_fold_pa, postList_fold_abund, XData_val, dfPi_val)
    gc()
  }
  
  if(!import_success) next
  
  # --- 6. Evaluate Hurdle Model Fit ---
  message("    Evaluating Hurdle Model Fit for ", label)
  MFEVAL <- evaluateModelFitRobust(hM = hM_full_abund, predY = predY)
  
  saveRDS(MFEVAL, file = mf_output_path)
  message("    Successfully saved Hurdle evaluation to: ", mf_output_path)
  
  rm(predY, MFEVAL)
  gc()
}

message("    End time: ", timestamp())