#' DESCRIPTION
#' S02b_fit_models_gpu_import_posterior.R: Reads and converts GPU-sampled posterior
#' chains back into fitted Hmsc model objects across all MCMC parameter configurations.
#'
#' METHODOLOGY
#' - Sources pipeline run configuration and posterior import utility functions.
#' - Iterates through each parameter grid row in the MCMC configuration.
#' - Reconstructs model run names and imports posterior samples into Hmsc model structures.
#'
#' INPUTS
#' - code/config_model.R: Central configuration parameters (model_id, mcmc configuration grid).
#' - code/_utilities_hmsc_gpu.R: GPU/HPC helper routines including `generate_run_name` and `import_posterior`.
#' - Posterior sample output files produced by GPU/Python sampling runs.
#'
#' OUTPUTS
#' - Reassembled fitted Hmsc model objects/files containing imported posterior distributions.

library(jsonify)
library(Hmsc)
library(tidyverse)

# --- 1. Load Configuration --- 
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")

# --- 2. Main Loop: Iterate Over the MCMC configuration ---
if (!run_config$model_mode %in% c("paired", "legacy")) {
  stop("model_mode must be 'paired' or 'legacy'.")
}
components <- if (run_config$model_mode == "paired") names(run_config$component_mcmc) else ""

for (component in components) {
  component_config <- run_config
  component_config$model_id <- paste0(run_config$model_id, component)
  if (run_config$model_mode == "paired") {
    component_config$mcmc <- run_config$component_mcmc[[component]]
  }
  message(paste0("--- Import for Model ID: ", component_config$model_id))

  for (i in 1:nrow(component_config$mcmc)) {
  #i <- 1
  run_name <- generate_run_name(component_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  mcmc_i <- component_config$mcmc[i, ]
  
  # ---B import posterior ---
  import_posterior(
    mcmc = mcmc_i,
    run_nm = run_name,
    config = component_config
  )
  }
}

