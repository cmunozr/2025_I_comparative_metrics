# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")

# --- 2. Configuration and Setup ---
run_name <- generate_run_name(run_config)
paths <- list(
  local_dir = getwd(),
  models_dir = file.path(getwd(), "models"),
  unfitted_models_file = file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
)
dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Load the Unfitted Model ---
load(file = paths$unfitted_models_file)
unfitted_model <- models[[run_config$model_id]]

# --- 4. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 5. Initialize Aggregator Lists ---
# These will store ALL commands from ALL MCMC loops
all_commands_aggregated <- list()
gpu_commands_aggregated <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)

# --- 6. Loop and Generate Commands ---
for(i in 1:nrow(mcmc_params)){
  
  mcmc_params_i <- mcmc_params[i, ]
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  base_model_name <- run_name[i] 
  
  run_specific_dir_local <- file.path(paths$models_dir, base_model_name)
  run_specific_dir_server <- file.path(run_config$server$server_models_dir, base_model_name)
  output_rds_path_local <- file.path(run_specific_dir_local, paste0(base_model_name, ".rds"))
  
  # Step A: Prepare the R-side Hmsc object
  prepared_model <- prepare_hpc_model(unfitted_model, mcmc_params_i)
  
  # Step B: Save the prepared model as JSON-RDS
  model_saved <- save_prepared_model(prepared_model, output_rds_path_local)
  
  # Step C: Generate commands and add to aggregators
  if (model_saved) {
    new_commands <- generate_commands(
      base_model_name = base_model_name, # Use the unique name
      run_specific_dir_server = run_specific_dir_server,
      mcmc_params = mcmc_params_i,
      run_config = run_config
    )
    
    # Aggregate Mode 1 and 2 commands
    all_commands_aggregated <- c(all_commands_aggregated, new_commands$all_commands)
    
    # Aggregate Mode 3 commands
    for (gpu_i in 1:run_config$gpu$n_gpus_available) {
      gpu_commands_aggregated[[gpu_i]] <- rbind(
        gpu_commands_aggregated[[gpu_i]], 
        new_commands$gpu_commands[[gpu_i]]
      )
    }
  }
}

# --- 7. Write .sh command scripts for python ---

write_commands_scripts(
  execution_mode = run_config$gpu$execution_mode, 
  txt_commands = unlist(all_commands_aggregated),
  gpu_commands = gpu_commands_aggregated,
  output_script_dir = paths$models_dir, 
  base_model_name = run_config$model_id,
  run_config = run_config
)

