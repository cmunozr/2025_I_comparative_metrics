# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(110724)

# --- 2. Configuration and Setup ---
run_name <- generate_run_name(run_config)
paths <- list(
  local_dir = getwd(),
  models_dir = file.path(getwd(), "models"),
  unfitted_models_file = file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
)
dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 4. Initialize Aggregator Lists ---
all_commands_aggregated <- list()
gpu_commands_aggregated <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)

# --- 5. Loop and Generate partition Commands ---
for(i in 1:nrow(mcmc_params)){
  # i <- 1
  mcmc_params_i <- mcmc_params[i, ]
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  base_model_name <- run_name[i] 
  
  # --- Load the fitted Model ---
  fitted_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  hM <- readRDS(fitted_model_path)

  # --- Make partition ----
  partition <- createPartition(hM, nfolds = run_config$cv$k)
  
  parts <- sort(unique(partition))
  
  # --- Create CV models ----
  hM_cv <- lapply(X = parts, FUN = function(X) set_cv_training_model(k = X, hM = hM, partition = partition))
  
  run_specific_dir_local <- file.path(paths$models_dir, base_model_name)
  run_specific_dir_server <- file.path(run_config$server$server_models_dir, base_model_name)
  
  cv_dir <- file.path(run_specific_dir_local, "cv")
  dir.create(cv_dir, showWarnings = F)
  
  cv_name <- paste0(base_model_name, "_cv_", parts)
  output_rds_path_local <- file.path(cv_dir, paste0(cv_name, ".rds"))
  
  saveRDS(partition, file = file.path(cv_dir, paste0("partition.rds")))
  
  for(p in 1:length(parts)){
    # p <- 1
    hm_cv_p <- hM_cv[[p]]
    
    # Step A: Prepare the R-side Hmsc object
    prepared_model <- prepare_hpc_model(hm_cv_p, mcmc_params_i)
    
    # Step B: Save the prepared model as JSON-RDS
    model_saved <- save_prepared_model(prepared_model, output_rds_path_local[p])
    
    # Step C: Generate commands and add to aggregators
    if (model_saved) {
      commands <- generate_commands(
        base_model_name = cv_name[p], # Use the unique name
        run_specific_dir_server = run_specific_dir_server,
        mcmc_params = mcmc_params_i,
        run_config = run_config
      )
      
      # Aggregate Mode 1 and 2 commands
      all_commands_aggregated <- c(all_commands_aggregated, commands$all_commands)
      
      # Aggregate Mode 3 commands
      for (gpu_i in 1:run_config$gpu$n_gpus_available) {
        gpu_commands_aggregated[[gpu_i]] <- rbind(
          gpu_commands_aggregated[[gpu_i]], 
          commands$gpu_commands[[gpu_i]]
        )
      }
    }
  }
  
}

# --- 7. Write .sh command scripts for python ---

write_commands_scripts(
  execution_mode = run_config$gpu$execution_mode, 
  txt_commands = unlist(all_commands_aggregated),
  gpu_commands = gpu_commands_aggregated,
  output_script_dir = paths$models_dir, 
  base_model_name = paste0(run_config$model_id, "_cv"),
  run_config = run_config
)

