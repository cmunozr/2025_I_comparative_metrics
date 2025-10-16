# This is the central configuration file for a single experimental run.

run_config <- list(
  
  # A unique ID for the model's structure (data, formulas, random effects).
  # This should match an entry in experiments_log/model_definition_log.csv
  model_id = "fbs_M002.1",
  
  # MCMC sampling parameters
  mcmc = data.frame(
    samples = 1000,
    thin = 1000,
    n_chains = 4,
    transient_proportion = 0.5, # following standard method on Hmsc course
    adapt_nf_proportion = 0.4 # following standard method on Hmsc course
  ),
  
  # GPU parameters
  gpu = list(
    execution_mode = 3,
    n_gpus_available = 2
  ),
  
  # Server environment (for execution_mode = 3) 
  server = list(
    python_env = "/home/avesta/munozcs/Documents/hmsc-venv",
    nohup_dir = "/home/avesta/munozcs/Documents/nohup_chains",
    # Define the base server path for model files
    server_models_dir = "models" 
  )
)

# Function to create the run name
generate_run_name <- function(config) {
  paste0(
    config$model_id,
    "_thin_", config$mcmc$thin,
    "_samples_", config$mcmc$samples,
    "_chains_", config$mcmc$n_chains
  )
}
