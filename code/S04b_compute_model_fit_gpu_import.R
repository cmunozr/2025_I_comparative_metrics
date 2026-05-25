# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# Define all strategies required for diagnosis (matches S04a)
validation_strategies <- c("metso_holdout", "route_blocked_cv", "random_cv")

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 4. Main Loop: Iterate Over MCMC config and Strategies ---
message(paste0("--- Import for Model ID: ", run_config$model_id))

for(i in 1:nrow(mcmc_params)){
  # Note: Fixed the indexing from S04a so base_model_name is correctly assigned
  run_name <- generate_run_name(run_config)[i]
  base_model_name <- run_name 
  mcmc_params_i <- mcmc_params[i, ]
  
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  # Load the full fitted model (used as the template for out-of-fold evaluation)
  fitted_model_path <- file.path("models", base_model_name, paste0("fitted_", base_model_name, ".rds"))
  if(!file.exists(fitted_model_path)) {
     message("  Full fitted model not found. Skipping to next MCMC config.")
     next
  }
  hM_full <- readRDS(fitted_model_path)
  
  for(strategy in validation_strategies) {
    message(paste0("\n  Processing strategy: ", strategy))
    
    # 4.1 Assign the correct label
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
    
    # 4.3 Import Posteriors for all folds
    models_list <- list()
    import_success <- TRUE
    
    for(p in 1:length(parts)){
      message("    Importing fold: ", parts[p])
      
      # import_posterior handles the JSON to RDS conversion
      res_path <- import_posterior(
        mcmc = mcmc_params_i, 
        config = run_config, 
        run_nm = base_model_name, 
        partition_number = parts[p], 
        label = label
      )
      
      if(is.logical(res_path) && res_path == FALSE) {
        warning("    Failed to import fold ", parts[p], ". Check Python logs. Skipping strategy.")
        import_success <- FALSE
        break
      } else {
         models_list[[p]] <- readRDS(res_path)
      }
    }
    
    if(!import_success) next
    
    # 4.4 Compute Predictions & Evaluate Model Fit
    message("    Calculating cross-validation predictions for ", label)
    gc() # Clean memory
    
    # Compute out-of-fold predicted values using the list of fitted fold-models
    predY <- Hmsc::computePredictedValues(hM_full, partition = partition, m = models_list, expected = TRUE)
    
    message("    Evaluating Model Fit for ", label)
    MFEVAL <- Hmsc::evaluateModelFit(hM = hM_full, predY = predY)
    
    # Save the Evaluation Object
    mf_output_path <- file.path(eval_dir, paste0("MF_", base_model_name, "_", label, ".rds"))
    saveRDS(MFEVAL, file = mf_output_path)
    message("    Successfully saved evaluation to: ", mf_output_path)
    
    # Clean memory before next strategy
    rm(models_list, predY, MFEVAL)
    gc()
  }
}