library(jsonify)
library(Hmsc)
library(tidyverse)

# 1. Load Configuration 
source("code/config_model.R") 

# 2. Main Loop: Iterate Over the MCMC configuration
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  #i <- 1
  
  # --- A. Setup for the Current Run ---
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  fitted_filename <- file.path("models", run_name, paste0("fitted_", run_name, ".RData"))
  
  # has this run already been imported?
  if (file.exists(fitted_filename)) {
    message("Output file already exists. Skipping.")
    next
  }
  
  # --- B. Find and Verify Posterior Chain Files ---
  run_specific_dir <- file.path("models", run_name)
  if (!dir.exists(run_specific_dir)) {
    warning(paste("Directory not found for", run_name, ". Skipping."))
    next
  }
  
  current_params <- run_config$mcmc[i, ]
  chains <- current_params$n_chains
  chain_indices <- 0:(chains - 1)
  expected_filenames <- paste0(run_name, "_post_chain", sprintf("%.2d", chain_indices), "_file.rds")
  chain_paths <- file.path(run_specific_dir, expected_filenames)
  
  existing_chains <- chain_paths[file.exists(chain_paths)]
  n_existing_chains <- length(existing_chains)
  
  if (length(existing_chains) < chains) {
    warning(paste0("Found ", length(existing_chains), " of ", chains, " expected chains for ", run_name, ". Take in mind!"))
  }
  
  # --- C. Import and Assemble the Model ---
  tryCatch({
    
    # Load the original unfitted model structure
    unfitted_rdata_path <- file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
    load(unfitted_rdata_path)
    unfitted_model <- models[[1]]
    
    post_list <- list()
    for (chain_path in existing_chains) {
      post_list[[length(post_list) + 1]] <- from_json(readRDS(chain_path)[[1]])[[1]]
    }
    
    samples <- current_params$samples
    thin <- current_params$thin
    transient <- ceiling(current_params$transient_proportion * samples * thin)
    
    fitted_model <- importPosteriorFromHPC(
      m = unfitted_model,
      postList = post_list,
      nSamples = samples,
      thin = thin,
      transient = transient,
      alignPost = TRUE
    )
    
    for(cInd in 1:n_existing_chains){
      if(unfitted_model$nr==1){
        for(i in 1:samples){
          fitted_model$postList[[cInd]][[i]]$Alpha = list(fitted_model$postList[[cInd]][[i]]$Alpha[1,])
        }
      }
      if(is.matrix(fitted_model$postList[[cInd]][[1]]$Alpha)){
        for(i in 1:samples){
          x = fitted_model$postList[[cInd]][[i]]$Alpha
          fitted_model$postList[[cInd]][[i]]$Alpha = lapply(seq_len(nrow(x)), function(i) x[i,])
        }
      }
    }
    
    # --- D. Save the Final, Individual Fitted Model Object ---
    save(fitted_model, file = fitted_filename)
    message("Successfully assembled and saved fitted model to: ", fitted_filename)
    
  }, error = function(e) {
    warning(paste("An error occurred while importing", run_name, ":", e$message))
  })
}

