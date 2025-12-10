cat("-----------START--------------\n")
library(Hmsc)
library(sf)
library(here)
library(doParallel)
library(foreach)
library(tidyverse)


source(file.path("code", "config_model.R"))

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# Turn off S2 geometry to avoid spherical geometry errors with EPSG:4326
sf::sf_use_s2(FALSE)
set.seed(11072024)

# 1. SETUP & DATA LOADING
cat(sprintf("[%s] Starting setup...\n", Sys.time()))
sufix <- "metso"

# Create output directory
pred_dir <- file.path("results", "predictions", sufix)
if(!dir.exists(pred_dir)) dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)


data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)
XData <- XData_list$XData |> 
  as.data.frame()

if(sufix == "metso"){
  # Load Covariates
  samp <- sample(seq_len(length(XData_list$polygon_id)), size = 2000)
  XData <- XData[samp, ]
  samp_standid <- XData_list$polygon_id[samp]
  
  # Load Spatial Data
  # MAINTAINING EPSG:4326 to match training data structure
  sp_df <- read_sf(here("data", "metso", "treatment_control_stand_filtered.gpkg")) |>
    st_centroid() |>
    st_transform("EPSG:4326") |>
    dplyr::filter(metso == 1, standid %in% XData_list$polygon_id) |>
    dplyr::distinct(standid, .keep_all = TRUE) |>
    dplyr::arrange(factor(standid, levels = XData_list$polygon_id))
  sp_df <- sp_df[samp,]
}else{
  sp_df <- read_sf(file.path("data", "metso", "donut_matches_sp.gpkg")) |> 
    st_centroid() |>
    st_transform("EPSG:4326") |>
    dplyr::filter(standid %in% XData_list$polygon_id) |>
    dplyr::distinct(standid, .keep_all = TRUE)
}

if(sufix == "metso") saveRDS(samp_standid, file.path(pred_dir, "pred_ids.rds"))
if(sufix == "control") saveRDS(sp_df$standid, file.path(pred_dir, "pred_ids.rds"))

# Load Model
run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)

# Extract Spatial Random Level Name
# This ensures the list passed to prepareGradient is named correctly (e.g., "route" or "Route")
spatial_level_name <- hM$rLNames[1] 
cat(sprintf("Spatial random level identified as: %s\n", spatial_level_name))

# 2. PARALLEL CONFIGURATION

# Configuration for batching
total_rows <- nrow(XData)
batch_size <- 100 #200
num_batches <- ceiling(total_rows / batch_size)

# Parallel Workers Setup
n_cores <- 20 #40 
cat(sprintf("Initializing parallel cluster with %d cores...\n", n_cores))

cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)

cat(sprintf("[%s] Starting PARALLEL prediction for %d stands in %d batches.\n", Sys.time(), total_rows, num_batches))

# 3. PARALLEL BATCH LOOP

parallel_time <- system.time({
  
  foreach(i = 1:num_batches, 
          .packages = c("Hmsc", "sf"), 
          .export = c("XData", "sp_df", "hM", "run_config", "spatial_level_name", "pred_dir", "batch_size", "total_rows")) %dopar% {
            
            # A. Define indices for this batch
            start_idx <- (i - 1) * batch_size + 1
            end_idx <- min(i * batch_size, total_rows)
            indices <- start_idx:end_idx
            
            cat(sprintf("[%s] Worker %d: Processing Batch %d/%d (Rows %d-%d)...\n", 
                        Sys.time(), Sys.getpid(), i, num_batches, start_idx, end_idx))
            
            # Define filename using model_id from your config
            batch_filename <- paste0("predY_", run_config$model_id, "_batch_", i, ".rds")
            batch_file_path <- file.path(pred_dir, batch_filename)
            
            # B. Check if file exists (Skip if done to allow restarting)
            if (!file.exists(batch_file_path)) {
              
              # C. Prepare Data Slices
              # Slice covariates
              XData_sub <- XData[indices, ]
              
              # Slice coordinates
              coords_sub <- st_coordinates(sp_df)[indices, ]
              
              # Create Named List for Gradient (Crucial for prepareGradient)
              sDataNew_sub <- list()
              sDataNew_sub[[spatial_level_name]] <- coords_sub
              
              # D. Prediction Block with Error Handling
              tryCatch({
                # 1. Prepare Gradient object
                Gradient <- prepareGradient(hM, XDataNew = XData_sub, sDataNew = sDataNew_sub)
                
                # 2. Predict 
                # expected = FALSE gives Integers (Abundance/Counts)
                predY <- predict(hM, Gradient = Gradient, expected = FALSE) 
                
                # 3. Save Result
                saveRDS(predY, file = batch_file_path)
                
                # Explicit cleanup to keep worker RAM usage stable
                rm(predY, Gradient, XData_sub, sDataNew_sub)
                gc()
                
                cat(sprintf("[%s] Worker %d: Batch %d COMPLETED.\n", Sys.time(), Sys.getpid(), i))
                
              }, error = function(e) {
                err_msg <- as.character(e)
                cat(sprintf("[%s] Worker %d: ERROR in Batch %d: %s\n", Sys.time(), Sys.getpid(), i, err_msg))
                writeLines(err_msg, file.path(pred_dir, paste0("ERROR_batch_", i, ".txt")))
              })
            }else{
              cat(sprintf("[%s] Worker %d: Batch %d SKIPPED (Exists).\n", Sys.time(), Sys.getpid(), i))
            }
            
            # Return NULL to minimize traffic back to the master process
            NULL 
          }
})

# Stop the cluster
stopCluster(cl)

cat("-----------DONE---------------\n")
print(parallel_time)