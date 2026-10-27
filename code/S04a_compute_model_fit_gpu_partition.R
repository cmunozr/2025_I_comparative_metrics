#' DESCRIPTION
#' S04a_compute_model_fit_gpu_partition.R: Generates spatial, holdout, or blocked
#' cross-validation partitions for fitted Hmsc models and builds GPU execution
#' scripts for fitting partitioned training models.
#'
#' METHODOLOGY
#' - Loads model run configurations and GPU computing utilities.
#' - Reads full fitted Hmsc models across MCMC grid configurations.
#' - Constructs partition based on configured validation strategies
#'   (e.g., METSO holdout, regional North-South split, route-blocked CV, or random CV).
#' - Creates training model subsets using set_training_model and serializes partition RDS files.
#' - Prepares models for HPC/GPU execution via JSON-RDS and generates dispatch shell scripts.
#'
#' INPUTS
#' - code/config_model.R: Model configuration parameters, MCMC grid, CV settings, and server directories.
#' - code/_utilities_hmsc_gpu.R: GPU helper functions (prepare_hpc_model, save_prepared_model, generate_commands, write_commands_scripts).
#' - models/{run_name}/fitted_{run_name}.rds: Fitted Hmsc model objects.
#' - {run_config$test$test_dir}: CSV containing holdout route flags (vakio, is_metso) for METSO holdout.
#' - data/fbs/vakiolinja/Vakiolinjat_routes_2393.gpkg: Spatial vector routes for geographic partitioning.
#' - data/ELY/ELY_Regions.shp: Administrative/regional boundaries for North-South partitioning.
#'
#' OUTPUTS
#' - models/{run_name}/{label}/partition.rds: Partition vectors indicating fold/holdout assignments.
#' - models/{run_name}/{label}/unfitted_{run_name}_{label}_{part}.rds: Unfitted Hmsc training models.
#' - models/{run_name}/{label}/{run_name}_{label}_{part}.rds: Prepared JSON-RDS models for HPC execution.
#' - models/{model_id}_diagnosis*.sh: Shell scripts for GPU job dispatch across partition folds.

# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(dplyr)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(11072024)

# Define all strategies required for diagnosis
validation_strategies <- c("metso_holdout", "north_south") # c("metso_holdout", "route_blocked_cv", "random_cv") # c("route_blocked_cv") #

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# --- 4. Initialize Aggregators ---
all_commands_aggregated <- list()
gpu_commands_aggregated <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)
validation_metadata <- list() # To track generated setups

# --- 5. Loop and Generate partition Commands ---
for(i in 1:nrow(mcmc_params)){
  # i <- 1
  run_name <- generate_run_name(run_config)[i]
  mcmc_params_i <- mcmc_params[i, ]
  message(paste0("\nProcessing config ", i, ": thin = ", mcmc_params_i$thin, ", samples = ", mcmc_params_i$samples))
  
  base_model_name <- run_name[i] 
  
  fitted_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  m <- readRDS(fitted_model_path)
  
  run_specific_dir_local <- file.path(models_dir, base_model_name)
  run_specific_dir_server <- file.path(run_config$server$server_models_dir, base_model_name)
  
  for(strategy in validation_strategies) {
    message(paste0("  Generating partition for strategy: ", strategy))
    
    if(strategy == "metso_holdout"){
      label <- "ho_metso"
      m$studyDesign$vakio <- m$studyDesign$vakio |> 
        as.character() |> 
        as.numeric()
      ho_routes <- read.csv(run_config$test$test_dir) |> 
        dplyr::select(vakio, is_metso) |> 
        dplyr::distinct(vakio, .keep_all = TRUE)
      study_design_with_ho <- m$studyDesign |>
        dplyr::left_join(ho_routes, by = "vakio")
      
      partition <- ifelse(
        is.na(study_design_with_ho$is_metso) | study_design_with_ho$is_metso == 0,
        1, 2
      )
      parts <- 1
      hM_ <- list(set_training_model(k = 2, hM = m, partition = partition))
      
    } else if (strategy == "route_blocked_cv") {
      label <- "cv_route"
      partition <- createPartition(m, nfolds = run_config$cv$k, column = "vakio")
      parts <- sort(unique(partition))
      hM_ <- lapply(X = parts, FUN = function(X) set_training_model(k = X, hM = m, partition = partition))
      
    } else if (strategy == "random_cv") {
      label <- "cv_random"
      partition <- createPartition(m, nfolds = run_config$cv$k)
      parts <- sort(unique(partition))
      hM_ <- lapply(X = parts, FUN = function(X) set_training_model(k = X, hM = m, partition = partition))
      
    } else if ( strategy == "north_south"){
      label <- "north_south"
      library(sf)
      geo_data_fbs <- read_sf("data/fbs/vakiolinja/Vakiolinjat_routes_2393.gpkg") |> 
        mutate(vakio = as.numeric(stringr::str_split(name,pattern = ",", simplify = T)[,1])) |> 
        select(vakio)
      geo_data_reg <- read_sf("data/ELY/ELY_Regions.shp") |> 
        st_transform(st_crs(geo_data_fbs))|> 
        mutate(north = if_else(ely == 15, 1, 0)) |> 
        select(north) |> 
        group_by(north) |> 
        summarise()
      
      fbs_with_regions <- st_join(geo_data_fbs, geo_data_reg) |> 
        st_drop_geometry()
      
      m$studyDesign$vakio <- m$studyDesign$vakio |> 
        as.character() |> 
        as.numeric()
      
      study_design_with_north <- m$studyDesign |>
        dplyr::left_join(fbs_with_regions, by = "vakio")
      
      partition <- ifelse(
        study_design_with_north$north == 0,
        1, 2
      )
      parts <- 1
      hM_ <- list(set_training_model(k = 2, hM = m, partition = partition))
      
    }
    
    # Save unfitted models
    dir.create(file.path("models", run_name, label), showWarnings = F, recursive = T)
    for(idx in 1:length(parts)){
      saveRDS(hM_[[idx]], file.path("models", run_name, label, paste0("unfitted_", run_name, "_", label, "_", parts[idx], ".rds")))
    }
    
    eval_dir <- file.path(run_specific_dir_local, label)
    dir.create(eval_dir, showWarnings = F)
    
    # Save the partition vector once for the label
    saveRDS(partition, file = file.path(eval_dir, "partition.rds"))
    
    # Record metadata
    validation_metadata[[length(validation_metadata) + 1]] <- data.frame(
      model_id = run_config$model_id,
      run_name = base_model_name,
      mcmc_row = i,
      strategy = strategy,
      label = label,
      nfolds = length(parts)
    )
    
    for(p in 1:length(parts)){
      # Define the unique name and path for the specific fold (p)
      eval_name_p <- paste0(base_model_name, "_", label, "_", parts[p])
      output_rds_path_local_p <- file.path(eval_dir, paste0(eval_name_p, ".rds"))
      
      prepared_model <- prepare_hpc_model(hM_[[p]], run_config$cv$mcmc)
      
      # Use the dynamically generated path
      model_saved <- save_prepared_model(prepared_model, output_rds_path_local_p) 
      
      if (model_saved) {
        commands <- generate_commands(
          base_model_name = eval_name_p, 
          run_specific_dir_server = file.path(run_specific_dir_server, label),
          mcmc_params = run_config$cv$mcmc,
          run_config = run_config
        )
        all_commands_aggregated <- c(all_commands_aggregated, commands$all_commands)
        for (gpu_i in 1:run_config$gpu$n_gpus_available) {
          gpu_commands_aggregated[[gpu_i]] <- rbind(gpu_commands_aggregated[[gpu_i]], commands$gpu_commands[[gpu_i]])
        }
      }
    } 
  }
}

# --- 6. Write Command Scripts ---

write_commands_scripts(
  execution_mode = run_config$gpu$execution_mode, 
  txt_commands = unlist(all_commands_aggregated),
  gpu_commands = gpu_commands_aggregated,
  output_script_dir = models_dir, 
  base_model_name = paste0(run_config$model_id, "_diagnosis"),
  run_config = run_config
)
