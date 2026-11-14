# This is the central configuration file for a single experimental run.

run_config <- list(
  
  # A unique ID for the model's structure (data, formulas, random effects).
  # This should match an entry in experiments_log/model_definition_log.csv
  model_id = "fbs_M017",
  model_mode = "paired", # "paired" or "legacy"
  
  # MCMC sampling parameters
  mcmc = data.frame(
    samples = 1000,
    thin = 250,
    n_chains = 4,
    transient_proportion = 0.5, # following standard method on Hmsc course
    adapt_nf_proportion = 0.4 # following standard method on Hmsc course
  ),

  component_mcmc = list(
    PA = data.frame(
      samples = 1000,
      thin = 150,
      n_chains = 4,
      transient_proportion = 0.5,
      adapt_nf_proportion = 0.4
    ),
    aCp = data.frame(
      samples = 1000,
      thin = 250,
      n_chains = 4,
      transient_proportion = 0.5,
      adapt_nf_proportion = 0.4
    )
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
  ),
  
  # Cross-validation
  cv = list(
    k = 4,
    mcmc = data.frame(
      samples = 1000,
      thin = 150,
      n_chains = 4,
      transient_proportion = 0.5, # following standard method on Hmsc course
      adapt_nf_proportion = 0.4 # following standard method on Hmsc course
    )
  ),
  # test-set
  test = list(
    test_dir = file.path("data", "fbs", "route_distance_2km_metso_control.csv")
  ) 
)

selected_models <- list(
  PA = "fbs_M017PA_thin_150_samples_1000_chains_4",
  aCp = "fbs_M017aCp_thin_250_samples_1000_chains_4"
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
