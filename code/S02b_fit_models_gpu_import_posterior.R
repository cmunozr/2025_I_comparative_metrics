library(jsonify)
library(Hmsc)
library(tidyverse)

# --- 1. Load Configuration --- 
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")

# --- 2. Main Loop: Iterate Over the MCMC configuration ---
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  #i <- 1
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  mcmc_i <- run_config$mcmc[i, ]
  
  # ---B import posterior ---
  import_posterior(
    mcmc = mcmc_i,
    run_nm = run_name,
    config = run_config
  )
}

