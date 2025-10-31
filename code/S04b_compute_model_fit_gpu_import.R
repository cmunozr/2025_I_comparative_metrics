# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(110724)

# --- 2. Configuration and Setup ---
paths <- list(
  local_dir = getwd(),
  models_dir = file.path(getwd(), "models"),
  unfitted_models_file = file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
)
dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc

# 4. Main Loop: Iterate Over the MCMC configuration
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  
  fitted_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  hM <- readRDS(fitted_model_path)
  
  #i <- 1
  parts <- 1:run_config$cv$k
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  mcmc_i <- run_config$mcmc[i, ]
  
  # --- A. Import to unfitted cross-validated model and read fitted cross-validation ---
  fitted_cv_paths <- lapply(
    X = parts, 
    FUN = function(X){
      import_posterior(
        mcmc = mcmc_i,
        run_nm = run_name,
        config = run_config, 
        partition_number = X
      )
    }
  ) 
  
  fitted_cv_models <- fitted_cv_paths |> 
    unlist() |> 
    lapply(readRDS)
  
  # --- B. Call partition ---
  base_model_name <- run_name[i] 
  run_specific_dir_local <- file.path(paths$models_dir, base_model_name)
  partition <- file.path(run_specific_dir_local, "cv", "partition.rds") |> 
    readRDS()
  
  # --- C. Combine predictions
  
  idfold <- parts
  partition.sp = NULL
  Yc = NULL
  expected=TRUE
  
  for(p in parts){
    #p <- 1
    message("predictions for partition ", p)
    val <- partition != p
    m <- do.call(c.Hmsc, mods[which(idfold == p)])
    m <- alignPosterior(m)
    
    postList = poolMcmcChains(hM1$postList, start=start, thin = thin)
    
    
  }
  
  
}
