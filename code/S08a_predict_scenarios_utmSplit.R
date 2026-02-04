cat("-----------START--------------\n")
library(Hmsc)
library(sf)
library(here)
library(doParallel)
library(foreach)
library(tidyverse)
library(arrow)
library(sfarrow)

source(file.path("code", "config_model.R"))
modelid <- run_config$model_id

if(!(Sys.getenv("RSTUDIO") == "1")){
  # setwd(here::here()) 
}

# Turn off S2 geometry to avoid spherical errors
sf::sf_use_s2(FALSE)
set.seed(11072024)

# --- CONFIGURATION ---
sampling <- FALSE        # <--- Set to FALSE for full run
sufix <- "control"       # <--- "metso" or "control"
expected_val <- TRUE
n_cores <- 20          # Number of cores
batch_size <- 4000       # Rows per batch, 200 in case of test
sampling_size <- 4000    # in case of test

# 1. SETUP & DATA LOADING

cat(sprintf("[%s] Starting setup for [%s]...\n", Sys.time(), sufix))

pred_dir <- file.path("results", "predictions")
if(!dir.exists(pred_dir)) dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

expected_string <- if(expected_val) "_expected_true" else "_expected_false"

pred_id_file <- file.path(dirname(pred_dir), paste0(sufix,"_pred_ids.rds"))

# Load Spatial Data
sp_df <- read_sf(here("data", "metso", "treatment_control_stand_v2.gpkg"))

# Categorizing by UTM zones
# Load UTM and Create Lookup
utm_200 <- st_read(file.path("data", "utm35_zones", "TM35_karttalehtijako.gpkg"), layer = "utm200")
utm_lookup <- sp_df |> 
  st_transform(st_crs(utm_200)) |> 
  st_join(utm_200) |> 
  st_drop_geometry() |> 
  group_by(standid) |> 
  slice(1) |> 
  ungroup() |> 
  dplyr::select(standid, lehtitunnus)

# Load Matches
matches <- readRDS(file.path("data", "metso", "raw", "matched_pairs.rds")) |> 
  filter(standid_treated %in% sp_df$standid | standid_matched_control %in% sp_df$standid) |>
  # We join UTM here just for the sampling/filtering logic if needed, 
  # but we will rejoin to sp_df later for the final grouping.
  left_join(utm_lookup, by = c("standid_treated" = "standid")) |> 
  rename(UTM200_metso = lehtitunnus) |> 
  left_join(utm_lookup, by = c("standid_matched_control" = "standid")) |> 
  rename(UTM200_control = lehtitunnus) |>
  mutate(same_UTM = UTM200_metso == UTM200_control) |> 
  na.omit()

saveRDS(matches, file.path("data", "metso", "matched_pairs_utm.rds"))


# Load XData
data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)
XData <- as.data.frame(XData_list$XData)

# Sampling Logic (Test Mode)
if(sufix == "metso" & sampling ){
  samp <- sample(seq_len(nrow(matches)), size = sampling_size)
  matches <- matches[samp,]
} 

if(sufix == "control"){
  # If prediction IDs already exist for metso, use them to filter control
  metso_ids_path <-  paste0(dirname(pred_id_file), "/metso_pred_ids.rds")
  if(file.exists(metso_ids_path)){
    samp_id <- readRDS(metso_ids_path)  
    matches <- matches |> filter(standid_treated %in% samp_id)
  }else{
    stop("First run metso and then control")
  }
}

# Define target IDs based on sufix
if(sufix == "metso") {
  target_ids <- matches$standid_treated
  target_metso <- 1
} else {
  target_ids <- matches$standid_matched_control
  target_metso <- 0
}

# Prepare sp_df (Centroids & Transform)
sp_df <- sp_df |>
  dplyr::filter(metso == target_metso, standid %in% target_ids) |>
  dplyr::distinct(standid, .keep_all = TRUE) |> 
  st_centroid() |>
  st_transform("EPSG:4326")

## Align Data
valid_ids <- intersect(sp_df$standid, XData_list$polygon_id)
sp_df <- sp_df[sp_df$standid %in% valid_ids, ]

match_idx <- match(sp_df$standid, XData_list$polygon_id)
XData <- XData[match_idx, , drop = FALSE]

## Safety Check
stopifnot(all(XData_list$polygon_id[match_idx] == sp_df$standid))

# save geometries
sp_df |> 
  dplyr::select(standid) |> 
  st_write_parquet(file.path("results", paste0("sites_geometry_", sufix, ".parquet")))

## Save IDs used for this run
saveRDS(sp_df$standid, file = pred_id_file)

coords <- st_coordinates(sp_df) |> as.data.frame()

# Attach utm info to final sp_df
# We join the lookup again to ensure every row in sp_df has a UTM zone
sp_df <- left_join(sp_df, utm_lookup, by = "standid")

# Load Model
run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)
spatial_level_name <- hM$rLNames[1] 

# 2. CREATE TASK LIST (Organized by UTM)

cat(sprintf("[%s] Generating Task List organized by UTM Zones...\n", Sys.time()))

unique_zones <- unique(na.omit(sp_df$lehtitunnus))
tasks <- list()
counter <- 1

for(zone in unique_zones) {
  # Get row indices for this specific UTM zone
  # Note: These indices correspond to the rows in XData/coords/sp_df
  zone_indices <- which(sp_df$lehtitunnus == zone)
  n_in_zone <- length(zone_indices)
  
  # Calculate how many batches are needed for this zone
  n_batches_zone <- ceiling(n_in_zone / batch_size)
  
  for(b in 1:n_batches_zone) {
    # Determine start/end within the zone_indices vector
    start_offset <- (b - 1) * batch_size + 1
    end_offset <- min(b * batch_size, n_in_zone)
    
    # Extract the actual global row indices
    batch_indices <- zone_indices[start_offset:end_offset]
    
    # Create Task Object
    tasks[[counter]] <- list(
      task_id = counter,
      utm_zone = zone,
      batch_num = b,            # Batch number WITHIN this zone
      total_batches_zone = n_batches_zone,
      indices = batch_indices,
      XData_batch = XData[batch_indices, ],
      coords_batch = coords[batch_indices, ],
      sp_df_batch = sp_df$standid[batch_indices]
    )
    
    counter <- counter + 1
  }
}

stopifnot(length(tasks) > 0)

cat(sprintf("Created %d tasks across %d UTM zones.\n", length(tasks), length(unique_zones)))

# 3. PARALLEL EXECUTION

cat(sprintf("Initializing parallel cluster with %d cores...\n", n_cores))
cl <- makeCluster(n_cores, outfile = "") # outfile="" prints worker output to console
registerDoParallel(cl)


parallel_time <- system.time({
  
  foreach(task = tasks, 
          .packages = c("Hmsc", "sf", "arrow", "dplyr", "tidyr", "tibble"), 
          .export = c("hM", "run_config", "modelid",
                      "spatial_level_name", "pred_dir",
                      "expected_val", "sufix", "expected_string")) %dopar% {
            
            # Unpack task info
            idx <- task$indices
            zone <- task$utm_zone
            b_num <- task$batch_num
            
            cat(sprintf("[%s] Worker %d: Processing UTM %s | Batch %d/%d (%d sites)...\n", 
                          Sys.time(), Sys.getpid(), zone, b_num, task$total_batches_zone, length(idx)))
              
            # Slice Data
            XData_sub <- task$XData_batch
            coords_sub <- task$coords_batch
              
            # Gradient List
            sDataNew_sub <- list()
            sDataNew_sub[[spatial_level_name]] <- coords_sub
            
            tryCatch({
              # Predict
              Gradient <- prepareGradient(hM, XDataNew = XData_sub, sDataNew = sDataNew_sub)
              predY <- predict(hM, Gradient = Gradient, expected = expected_val, predictEtaMean = TRUE) |> 
                simplify2array()
              
              # predY <- readRDS("results/predY_fbs_M008_batch_1.rds") |> simplify2array()
              
              batch_ids <- task$sp_df_batch
                
              # We need to convert the 3D array to a 2D Data Frame.
              # Strategy: Columns = standid, iteration, species_1, species_2.
              
              # A. Define Dimensions
              n_sites_batch <- dim(predY)[1]
              n_spp <- dim(predY)[2]        
              n_post <- dim(predY)[3]
              
              # B. Rotate the Array
              # Current: [Site, Species, Posterior]
              # Target:  [Site, Posterior, Species]
              predY <- aperm(predY, c(1, 3, 2))
              
              # C. Squash into Matrix. 
              # Combine Dim 1 (Site) and Dim 2 (Iter) into the rows. Dim 3 (Species) becomes columns.
              dim(predY) <- c(n_sites_batch * n_post, n_spp)
              
              # D. Convert to Data Frame (Lightweight)
              predY <- as.data.frame(predY)
              
              # Fix column names (Species names)
              colnames(predY) <- colnames(hM$Y)
              
              # E. Add ID Columns
              # Construct the IDs to match the order of the wide matrix
              # Row 1: Site 1, Iter 1
              # Row 2: Site 2, Iter 1
              # ...
              # Row 4000: Site 4000, Iter 1
              # Row 4001: Site 1, Iter 2
              
              site_ids_rep <- rep(batch_ids, times = n_post)
              post_ids_rep <- rep(1:n_post, each = n_sites_batch)
              
              # Attach IDs
              predY <- predY |>
                mutate(
                  standid = site_ids_rep,
                  posterior = post_ids_rep
                ) |>
                relocate(standid, posterior, year) 
              
              # --- MANUAL DATASET CREATION ---
              
              # Build the Hive Partition
              part_path <- file.path(
                pred_dir,
                paste0("modelid=", modelid),
                paste0("expected_type=", expected_string),
                paste0("year=2021"),
                paste0("utm_zone=", zone),
                paste0("scenario=", sufix)
              )
              
              # 2. Create the directory safely
              # recursive = TRUE ensures all parent folders are created
              if(!dir.exists(part_path)) {
                dir.create(part_path, recursive = TRUE, showWarnings = FALSE)
              }
              
              # Create unique filename
              file_name <- paste0("part_batch", b_num, ".parquet")
              full_file_path <- file.path(part_path, file_name)
              
              # Write the Parquet file
              arrow::write_parquet(predY, sink = full_file_path)
              
              # Cleanup massive objects immediately
              rm(site_ids_rep, post_ids_rep)
                
              }, error = function(e) {
                err_msg <- as.character(e)
                cat(sprintf("ERROR in UTM %s Batch %d: %s\n", zone, b_num, err_msg))
                writeLines(err_msg, file.path(pred_dir, paste0("ERROR_UTM_", sufix, zone, "_batch_", b_num, ".txt")))
              })
            
            NULL
          }
})

stopCluster(cl)
cat("-----------DONE---------------\n")
print(parallel_time)