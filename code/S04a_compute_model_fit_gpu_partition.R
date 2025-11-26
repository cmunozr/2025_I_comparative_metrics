# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# This is a Crossvalidation model fit or a spatial hold-out
do_spatial_holdout <- TRUE

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 4. Initialize Aggregator Lists ---
all_commands_aggregated <- list()
gpu_commands_aggregated <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)

# --- 5. Loop and Generate partition Commands ---
for(i in 1:nrow(mcmc_params)){
  # i <- 1
  run_name <- generate_run_name(run_config)[i]
  mcmc_params_i <- mcmc_params[i, ]
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  base_model_name <- run_name[i] 
  
  # --- Load the fitted Model ---
  fitted_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  m <- readRDS(fitted_model_path)
  
  run_specific_dir_local <- file.path(models_dir, base_model_name)
  run_specific_dir_server <- file.path(run_config$server$server_models_dir, base_model_name)
  
  if(do_spatial_holdout){
    
    label <- "ho"
    m$studyDesign$vakio <- as.numeric(m$studyDesign$vakio)
    
    # --- Make partition ----
    ho_routes <- read.csv(run_config$test$test_dir) |> 
      dplyr::select(vakio, is_metso) |> 
      dplyr::distinct(vakio, .keep_all = TRUE)
    
    study_design_with_ho <- m$studyDesign |>
      dplyr::left_join(ho_routes, by = "vakio")
    
    partition <- ifelse(
      !is.na(study_design_with_ho$is_metso) & study_design_with_ho$is_metso == 1,
      2, # Fold 2 = Test Set
      1  # Fold 1 = Training Set
    )
    
    parts <- 1
    
    # --- Create ho models ----
    hM_ <- list(set_training_model(k = 2, hM = m, partition = partition))
    dir.create(file.path("models", run_name, label), showWarnings = F, recursive = T)
    saveRDS(hM_, file.path("models", run_name, label, paste0("unfitted_", run_name, "_", label, "_", parts, ".rds")))
    
  }else{
    
    label <- "cv"
    
    # --- Make partition ----
    partition <- createPartition(m, nfolds = run_config$cv$k, column = "vakio")
    
    parts <- sort(unique(partition))
    
    # --- Create CV models ----
    hM_ <- lapply(X = parts, FUN = function(X){
      cv_t_model <- set_training_model(k = X, hM = m, partition = partition)
      dir.create(file.path("models", run_name, label), showWarnings = F, recursive = T)
      saveRDS(cv_t_model, file.path("models", run_name, label, paste0("unfitted_", run_name, "_", label, "_", X, ".rds")))
      return(cv_t_model)
    }) 
    
  }
  
  eval_dir <- file.path(run_specific_dir_local, label)
  dir.create(eval_dir, showWarnings = F)
  
  eval_name <- paste0(base_model_name, "_", label, "_", parts)
  output_rds_path_local <- file.path(eval_dir, paste0(eval_name, ".rds"))
  
  saveRDS(partition, file = file.path(eval_dir, paste0("partition.rds")))
  
  for(p in 1:length(parts)){
    # p <- 1
    hm_p <- hM_[[p]]
    
    # Step A: Prepare the R-side Hmsc object
    prepared_model <- prepare_hpc_model(hm_p, run_config$cv$mcmc_temp)
    
    # Step B: Save the prepared model as JSON-RDS
    model_saved <- save_prepared_model(prepared_model, output_rds_path_local[p])
    
    # Step C: Generate commands and add to aggregators
    if (model_saved) {
      commands <- generate_commands(
        base_model_name = eval_name[p], # Use the unique name
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

# --- 6. Write .sh command scripts for python ---

write_commands_scripts(
  execution_mode = run_config$gpu$execution_mode, 
  txt_commands = unlist(all_commands_aggregated),
  gpu_commands = gpu_commands_aggregated,
  output_script_dir = models_dir, 
  base_model_name = paste0(run_config$model_id, "_", label),
  run_config = run_config
)
