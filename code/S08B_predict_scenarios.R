#' DESCRIPTION:
#' This script generates batch-processed spatial predictions for HMSC models.
#' It splits treated (Metso) and control stands into UTM-based chunks to 
#' parallelize the heavy posterior sampling and saves the results into a 
#' partitioned Parquet data lake for easy comparison.
#'
#' METHODOLOGY:
#' - Parallelism: Uses 'doParallel' and 'foreach' for multi-core acceleration.
#' - Batching: Spatial data is segmented by UTM zones and sub-divided into 
#'   manageable row-counts (batch_size) to optimize memory overhead.
#' - Storage: Implements a Hive-style partitioned directory structure for 
#'   Parquet exports (modelid / expected_type / year / utm_zone / scenario).
#' - Geometry: Centroid-based prediction with S2 geometry disabled for 
#'   planar coordinate consistency.
#'
#' INPUTS:
#' - HMSC Model: Fitted model object (.rds) located via 'config_model.R'.
#' - Spatial Data: 'treatment_control_stand_v2.gpkg' and UTM zone lookups.
#' - Covariates: 'XData_hmsc_[sufix]_[model_id].rds'.
#' - Matched Pairs: 'matched_pairs.rds' for treated/control filtering.
#'
#' OUTPUTS:
#' - Predictions: Hive-partitioned Parquet files (one per batch).
#' - Geometries: 'sites_geometry_[sufix].parquet'.
#' - Metadata: 'matched_pairs_utm.rds' and prediction ID trackers.


library(Hmsc)
library(sf)
library(here)
library(doParallel)
library(foreach)
library(tidyverse)
library(arrow)
library(sfarrow)
library(abind)

source(file.path("code", "config_model.R"))
modelid <- run_config$model_id

if(!(Sys.getenv("RSTUDIO") == "1")){
    setwd(here::here()) 
}

# Turn off S2 geometry to avoid spherical errors
sf::sf_use_s2(FALSE)
set.seed(11072024)

# --- CONFIGURATION ---
sufix <- "metso" # or "control" (it should change to BAU in some point)
if(sufix == "metso"){
  val = 1
}else{
  val = 0
}

test <- TRUE             # <--- Set to FALSE for full run
expected_val <- TRUE

if(test){
  batch_size <- 100      
  sampling_size <- 1000      
  n_cores <- 2
}else{
  batch_size <- 4000
  n_cores <- 20
}

# alfa matrix
alfa_matrix <- readRDS(file.path("results", "alfa_matrix.rds"))

# Load presence/abscence Model
fitted_full_model_path_PA <- file.path("models/fbs_M016PA_thin_150_samples_1000_chains_4/fitted_fbs_M016PA_thin_150_samples_1000_chains_4.rds")
hM_PA <- readRDS(fitted_full_model_path_PA)
spatial_level_name <- hM_PA$rLNames[1] 

# Load abundance conditional on presence Model
fitted_full_model_path_aCp <- file.path("models/fbs_M016_thin_250_samples_1000_chains_4/fitted_fbs_M016_thin_250_samples_1000_chains_4.rds")
hM_aCp <- readRDS(fitted_full_model_path_aCp)

#-----------------------

# 1. SETUP & DATA LOADING

cat(sprintf("[%s] Starting setup \n", Sys.time()))

pred_dir <- file.path("results", "predictions")
if(!dir.exists(pred_dir)) dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

expected_string <- if(expected_val) "_expected_true" else "_expected_false"

pred_id_file <- file.path(dirname(pred_dir), paste0("pred_ids_", sufix, ".rds"))

# Load Spatial Data
sp_df <- read_sf(here("data", "metso", "treatment_control_stand_v4.gpkg")) |> 
  dplyr::filter(metso == val)

# create group Matches
matches <- sp_df |> 
  group_by(year, regional_group, treespecies) |> #ely_en
  mutate(group_number = cur_group_id()) |> 
  ungroup() |> 
  select(standid, metso, group_number)

# Load XData
data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)
XData <- as.data.frame(XData_list$XData)
XData_polygon <- XData_list$polygon_id

## Align Data
valid_ids <- intersect(sp_df$standid, XData_polygon)
sp_df <- sp_df[sp_df$standid %in% valid_ids, ]
matches <- matches[matches$standid %in% valid_ids, ]

match_idx <- match(sp_df$standid, XData_polygon)
XData <- XData[match_idx, , drop = FALSE]

## Safety Check
stopifnot(all(XData_polygon[match_idx] == sp_df$standid))

# transform and get geometries

sp_df <- sp_df |>
  dplyr::distinct(standid, .keep_all = TRUE) |> 
  st_centroid() |>
  st_transform("EPSG:4326")

coords <- st_coordinates(sp_df) |> 
  as.data.frame()

# Save matches
saveRDS(matches, file.path("results", paste0("matched_pairs_groups", "_", sufix, ".rds")))

# Save IDs used for this run
saveRDS(sp_df$standid, file = pred_id_file)

# save geometries
sp_df |> 
  dplyr::select(standid) |> 
  st_write_parquet(file.path("results", paste0("point_geometry_", sufix, ".parquet")))

# test Logic (Test Mode), select just some sites
if(test ){
  samp <- sample(seq_len(nrow(matches)), size = sampling_size)
  matches <- matches[samp, ]
  sp_df <- sp_df[samp, ]
  coords <- coords[samp, ]
  XData <- XData[samp, ]
} 


# complete XData, Define reference values from the training data hM_PA$XData
XData_survey <- hM_PA$XData
missing_cols <- setdiff(colnames(XData_survey), colnames(XData))

ref_values <- list()
for (col in missing_cols) {
  vals <- XData_survey[[col]]
  if (is.numeric(vals)) {
    # Mean of the training survey variable
    ref_values[[col]] <- mean(vals, na.rm = TRUE)
  } else if (is.factor(vals)) {
    # Most common level or default reference level
    ref_values[[col]] <- levels(vals)[1]
  } else {
    ref_values[[col]] <- unique(vals)[1]
  }
}

# Function to expand any stand-level batch (XData_sub) with the constant reference survey covariates
prepare_prediction_xdata <- function(stand_data, ref_vals, ref_template) {
  df_out <- stand_data
  for (nm in names(ref_vals)) {
    df_out[[nm]] <- ref_vals[[nm]]
  }
  # Ensure column order and factor levels match the training set exactly
  for (nm in colnames(ref_template)) {
    if (is.factor(ref_template[[nm]])) {
      df_out[[nm]] <- factor(df_out[[nm]], levels = levels(ref_template[[nm]]))
    }
  }
  return(df_out)
}

XData <- prepare_prediction_xdata(
  stand_data = XData, 
  ref_vals = ref_values, 
  ref_template = XData_survey
)

#-----------------------------

# 2. CREATE TASK LIST (Organized by UTM)

cat(sprintf("[%s] Generating Task List organized by groups...\n", Sys.time()))

unique_groups <- unique(na.omit(matches$group_number))
tasks <- list()
counter <- 1

for(group in unique_groups) {
  # Get row indices for this specific group
  # Note: These indices correspond to the rows in XData/coords/sp_df
  group_indices <- which(matches$group_number == group)
  n_in_group <- length(group_indices)
  
  # Calculate how many batches are needed for this group
  n_batches_group <- ceiling(n_in_group / batch_size)
  
  for(b in 1:n_batches_group) {
    # Determine start/end within the group_indices vector
    start_offset <- (b - 1) * batch_size + 1
    end_offset <- min(b * batch_size, n_in_group)
    
    # Extract the actual global row indices
    batch_indices <- group_indices[start_offset:end_offset]
    
    # Create Task Object
    tasks[[counter]] <- list(
      task_id = counter,
      group = group,
      year = sp_df$year[batch_indices],
      regional = sp_df$regional_group[batch_indices],
      #ely = sp_df$ely_en[batch_indices],
      trees = sp_df$treespecies[batch_indices],
      batch_num = b,            # Batch number WITHIN this group
      total_batches_group = n_batches_group,
      indices = batch_indices,
      XData_batch = XData[batch_indices, ],
      coords_batch = coords[batch_indices, ],
      sp_df_batch = sp_df$standid[batch_indices]
    )
    
    counter <- counter + 1
  }
}

stopifnot(length(tasks) > 0)

cat(sprintf("Created %d tasks across %d groups.\n", length(tasks), length(unique_groups)))

# 3. PARALLEL EXECUTION

start_time <- Sys.time()

cat(sprintf("Initializing parallel cluster with %d cores...\n", n_cores))
cl <- makeCluster(n_cores, outfile = "") # outfile="" prints worker output to console
registerDoParallel(cl)


foreach(task = tasks, 
        .packages = c("Hmsc", "sf", "arrow", "dplyr", "tidyr", "tibble"), 
        .export = c("hM_PA", "hM_aCp", "run_config", "modelid",
                    "spatial_level_name", "pred_dir",
                    "expected_val", "sufix", "expected_string",
                    "alfa_matrix")) %dopar% {
            
          # Unpack task info
          idx <- task$indices
          group <- task$group
          yr <- task$year
          reg <- task$regional
          #ely <-  task$ely
          tr <-  task$trees
          b_num <- task$batch_num
            
          cat(sprintf("[%s] Worker %d: Processing UTM %s | Batch %d/%d (%d sites)...\n", 
                        Sys.time(), Sys.getpid(), group, b_num, task$total_batches_group, length(idx)))
              
          # Slice Data
          XData_sub <- task$XData_batch
          coords_sub <- task$coords_batch
              
          # Gradient List
          sDataNew_sub <- list()
          sDataNew_sub[[spatial_level_name]] <- coords_sub
            
          tryCatch({
            # Predict
            Gradient <- prepareGradient(hM_PA, XDataNew = XData_sub, sDataNew = sDataNew_sub)
            predY_PA <- predict(hM_PA, Gradient = Gradient, expected = expected_val, predictEtaMean = TRUE) |> 
              simplify2array()
            
            Gradient <- prepareGradient(hM_aCp, XDataNew = XData_sub, sDataNew = sDataNew_sub)
            predY_aCP <- predict(hM_aCp, Gradient = Gradient, expected = expected_val, predictEtaMean = TRUE) |> 
              simplify2array()
            
            n_sites_batch   <- dim(predY_PA)[1]
            n_spp      <- dim(predY_PA)[2]
            n_post <- dim(predY_PA)[3]
            
            predY_hurdle <- array(
              0,
              dim = c(n_sites_batch, n_spp, n_post),
              dimnames = dimnames(predY_PA)
            )
            
            for (s in seq_len(n_samples)) {
              
              th_s <- alfa_matrix[, s]
              is_present <- sweep(predY_PA[, , s], 2, th_s, ">=")
              
              # Apply threshold: keep probabilities where TRUE, zero out where FALSE, alfa
              pa_thresh <- predY_PA[, , s] * is_present
              
              # Construct hurdle: multiply thresholded probability by conditional abundance
              predY_hurdle[, , s] <- pa_thresh * predY_aCP[, , s]         
              
            }
            
            # predY <- readRDS("results/predY_fbs_M008_batch_1.rds") |> simplify2array()
              
            batch_ids <- task$sp_df_batch
                
            # We need to convert the 3D array to a 2D Data Frame.
            # Strategy: Columns = standid, iteration, species_1, species_2.
              
            # B. Rotate the Array
            # Current: [Site, Species, Posterior]
            # Target:  [Site, Posterior, Species]
            predY <- aperm(predY_hurdle, c(1, 3, 2))
              
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
              paste0("year=", yr),
              paste0("reg=", reg),
              #paste0("ely=", ely),
              paste0("trees=", tr),
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
            
            # Cleanup 
            rm(site_ids_rep, post_ids_rep)
                
            }, error = function(e) {
              err_msg <- as.character(e)
              cat(sprintf("ERROR in group %s Batch %d: %s\n", group, b_num, err_msg))
              writeLines(err_msg, file.path(pred_dir, paste0("ERROR_UTM_", sufix, group, "_batch_", b_num, ".txt")))
            })
          
          NULL
        }

stopCluster(cl)

end_time <- Sys.time()
total_time <- end_time - start_time

cat(sprintf("[%s] Done. Total time: %s\n", 
            Sys.time(), 
            format(total_time, digits = 2)))