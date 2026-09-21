#' DESCRIPTION:
#' This script generates parallelized, batch-processed spatial hurdle predictions 
#' for HMSC models across forest stands (Metso/treated or BAU/control).
#' It segments stands into stratified groups (year, regional_group, treespecies) 
#' and processes them in manageable site batches across multiple CPU cores. 
#' The resulting posterior hurdle predictions are saved into a Hive-partitioned 
#' Parquet dataset.
#'
#' METHODOLOGY:
#' - Hurdle Construction: Independent predictions are generated for PA (hM_PA) and 
#'   aCp (hM_aCp). For each posterior sample, species PA probabilities are compared 
#'   against an alpha threshold matrix (e.g., Minimum Training Presence, MTP); values 
#'   below threshold are zeroed out, PA probabilities >= alfa are multiplied 
#'   by conditional abundances to produce the hurdle prediction.
#' - Stratification & Batching: Stands are grouped by (year, regional_group, treespecies) 
#'   and chunked by 'batch_size' to optimize memory during posterior sampling.
#' - Covariate Alignment: Stand-level covariates are expanded to match model training 
#'   specifications by imputing missing survey covariates with mean baseline values.
#' - Geometry: Stand polygons are reduced to centroids, projected to WGS 84 (EPSG:4326).
#' - Storage: Predictions are written as individual batch Parquet files in a Hive-style 
#'   partition hierarchy: modelid / expected_type / year / reg / trees / scenario.
#' - Parallelism: Multi-core cluster execution handled via 'doParallel' and 'foreach'.
#'
#' INPUTS:
#' - HMSC Models: Fitted model RDS files for PA ('hM_PA') and aCp ('hM_aCp').
#' - Model Config: 'code/config_model.R' defining 'model_id'.
#' - Threshold Matrix: 'results/alfa_matrix.rds' (species x posterior thresholds).
#' - Spatial Data: Stand polygons from 'data/metso/treatment_control_stand_v4.gpkg'.
#' - Covariates: Extracted predictor RDS at 'data/covariates/XData_hmsc_[sufix]_[model_id].rds'.
#'
#' OUTPUTS:
#' - Predictions: Hive-partitioned Parquet files ('part_batch[B].parquet').
#' - Geometries: Centroid point layer saved as 'results/point_geometry_[sufix].parquet'.
#' - Metadata: 'matched_pairs_groups_[sufix].rds' and sampled IDs at 'pred_ids_[sufix].rds'.
#' - Logs: 'ERROR_group_[...].txt' files on failed worker tasks..


library(Hmsc)
library(sf)
library(here)
library(doParallel)
library(foreach)
library(tidyverse)
library(arrow)
library(sfarrow)
library(abind)
library(stringr)

source(file.path("code", "config_model.R"))
source(file.path("code", "_utilities_transform_covariates.R"))
source(file.path("code", "_utilities_diversity_metrics_functions.R"))
modelid <- run_config$model_id

message("    Start time: ", Sys.time())

if(!(Sys.getenv("RSTUDIO") == "1")){
    setwd(here::here()) 
}

# Turn off S2 geometry to avoid spherical errors
sf::sf_use_s2(FALSE)
set.seed(11072024)

# --- CONFIGURATION ---
sufix <- "control" # or "control" (it should change to BAU in some point)
if(sufix == "metso"){
  val = 1
}else{
  val = 0
}

# Execution mode flag: set to FALSE Set to FALSE for full run
test <- FALSE 

# Expected abundance predictions vs stochastic realizations (draws from distribution)
expected_val <- TRUE

if(test){
  batch_size <- 50      
  n_groups <- 20
  min_stands <- 5
  max_stands <- 50
}else{
  batch_size <- 500
}

os <- Sys.info()['sysname']

if (os == "Windows") {
  
  wr_dir <- "results"
  worker_tmp <- tempdir()
  n_cores <- max(1, floor(parallel::detectCores() / 2))
  
} else if (os %in% c("Linux", "Darwin")) {
  
  n_cores <- min(32, floor(parallel::detectCores() / 2))
  
  node <- Sys.info()["nodename"]
  
  if (node == "lema") {
    wr_dir <- NULL
    worker_tmp <- NULL
  } else if(node == "gurnah"){
    wr_dir <- "/export/scratch/tmp/munozcs/results"
    worker_tmp <- "/export/scratch/tmp/munozcs/Rtmp"
  } else {
    wr_dir <- NULL # here::here("results")
    worker_tmp <- NULL # tempdir()
  }
}

# Set to FALSE for full runs
retry_failed_only <- TRUE

#-----------------------

# 1. DATA LOADING

# Load presence/abscence Model
fitted_full_model_path_PA <- file.path("models/fbs_M016PA_thin_150_samples_1000_chains_4/fitted_fbs_M016PA_thin_150_samples_1000_chains_4.rds")
hM_PA <- readRDS(fitted_full_model_path_PA)
spatial_level_name <- hM_PA$rLNames[1] 

# Load abundance conditional on presence Model
fitted_full_model_path_aCp <- file.path("models/fbs_M016_thin_250_samples_1000_chains_4/fitted_fbs_M016_thin_250_samples_1000_chains_4.rds")
hM_aCp <- readRDS(fitted_full_model_path_aCp)

if(is.null(wr_dir)){
  stop("Warning: Check for routes to the physical disk without NFS mounting")
}

if(!dir.exists(wr_dir)) {
  stop("First run S08a or copy the outputs elements created there to the expected working directory /results/alfa_matrix.rds")
}

# alfa matrix
alfa_matrix <- readRDS(file.path(wr_dir, "alfa_matrix.rds"))

# output folder
pred_dir <- file.path(wr_dir, "predictions")

if(!dir.exists(pred_dir)) dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

expected_string <- if(expected_val) "_expected_true" else "_expected_false"

pred_id_file <- file.path(dirname(pred_dir), paste0("pred_ids_", sufix, ".rds"))

# Load Spatial Data
sp_df <- read_sf(here("data", "metso", "treatment_control_stand_v4.gpkg")) |> 
  dplyr::filter(metso == val) |> 
  dplyr::mutate(ely_en = stringi::stri_trans_general(ely_en, "Latin-ASCII"),
                ely_en = stringr::str_to_snake(ely_en))

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

if (test) {
  test_meta <- get_or_create_test_groups(
    gpkg_path = here("data", "metso", "treatment_control_stand_v4.gpkg"),
    test_group_path = file.path(wr_dir, "test_group_stands.rds"),
    n_groups = n_groups,
    min_stands = min_stands,
    max_stands = max_stands,
    seed = 11072024
  )
  
  # Extract IDs relevant for current scenario (sufix: metso or control)
  target_test_ids <- test_meta |>
    dplyr::filter(metso == val) |>
    dplyr::pull(standid)
  
  valid_test_ids <- intersect(valid_ids, target_test_ids)
  sp_df   <- sp_df[sp_df$standid %in% valid_test_ids, ]
  matches <- matches[matches$standid %in% valid_test_ids, ]
}

match_idx <- match(sp_df$standid, XData_polygon)
XData <- XData[match_idx, , drop = FALSE]

## Safety Check
stopifnot(all(XData_polygon[match_idx] == sp_df$standid))

# transform and get geometries

sp_df <- sp_df |>
  st_centroid() |>
  st_transform("EPSG:4326")

coords <- st_coordinates(sp_df) |> 
  as.data.frame()

# Save matches
saveRDS(matches, file.path(wr_dir, paste0("matched_pairs_groups", "_", sufix, ".rds")))

# Save IDs used for this run
saveRDS(sp_df$standid, file = pred_id_file)

# save geometries
sp_df |> 
  dplyr::select(standid) |> 
  st_write_parquet(file.path(wr_dir, paste0("point_geometry_", sufix, ".parquet")))

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

XData <- prepare_prediction_xdata(
  stand_data = XData, 
  ref_vals = ref_values, 
  ref_template = XData_survey,
  hM = hM_PA
)


#-----------------------------

# 2. CREATE TASK LIST (Organized by group)

cat(sprintf("[%s] Generating Task List organized by groups...\n", Sys.time()))

unique_groups <- unique(na.omit(matches$group_number)) |> sort()
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

#---------------------------
# Retry!

if (retry_failed_only) {
  tasks <- filter_failed_tasks(
    tasks            = tasks,
    pred_dir         = pred_dir,
    modelid          = modelid,
    expected_string  = expected_string,
    sufix            = sufix,
    mode             = "missing_output", 
    clean_error_logs = TRUE              # Deletes old ERROR_*.txt for retried tasks
  )
}

# If all tasks already completed successfully, avoid spinning up the cluster
if (length(tasks) == 0) {
  stop("All batches have already been processed successfully. Exiting.")
}

#--------------------------

# 3. PARALLEL EXECUTION

# Set TMPDIR for the main R session and Arrow
Sys.setenv(TMPDIR = worker_tmp)

cl <- makeCluster(n_cores, outfile = "")
cat(sprintf("Initializing parallel cluster with %d cores...\n", n_cores))
clusterExport(cl, "worker_tmp")
clusterEvalQ(cl, {
  if (!dir.exists(worker_tmp)) dir.create(worker_tmp, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TMPDIR = worker_tmp)
})

distinct_partitions <- sp_df |>
  sf::st_drop_geometry() |>
  distinct(year, regional_group, treespecies) #ely_en

for(i in 1:nrow(distinct_partitions)) {
  p_path <- file.path(pred_dir,
                      paste0("modelid=", modelid),
                      paste0("expected_type=", expected_string),
                      paste0("year=", distinct_partitions$year[i]),
                      paste0("reg=", distinct_partitions$regional_group[i]),
                      #paste0("ely=", distinct_partitions$ely_en[i]),
                      paste0("trees=", distinct_partitions$treespecies[i]),
                      paste0("scenario=", sufix)
  )
  dir.create(p_path, recursive = TRUE, showWarnings = FALSE)
}

registerDoParallel(cl)

foreach(task = tasks, 
        .packages = c("Hmsc", "sf", "arrow", "dplyr", "tidyr", "tibble"), 
        .export = c("hM_PA", "hM_aCp", "run_config", "modelid",
                    "spatial_level_name", "pred_dir",
                    "expected_val", "sufix", "expected_string",
                    "alfa_matrix")) %dopar% {
          
          # task <- tasks[[3]]  
          # Unpack task info
          idx <- task$indices
          group <- task$group
          yr <- task$year[1]
          reg <- task$regional[1]
          #ely <-  task$ely[1]
          tr <-  task$trees[1]
          b_num <- task$batch_num[1]
            
          cat(sprintf("[%s] Worker %d: Processing group %s | Batch %d/%d (%d sites)...\n", 
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
            predY_PA <- predict(object = hM_PA, X = as.matrix(Gradient$XDataNew), ranLevels = Gradient$rLNew, 
                                studyDesign = Gradient$studyDesignNew, expected = expected_val, predictEtaMean = TRUE) |> 
              simplify2array()
            
            Gradient <- prepareGradient(hM_aCp, XDataNew = XData_sub, sDataNew = sDataNew_sub)
            predY_aCP <-  predict(object = hM_aCp, X = as.matrix(Gradient$XDataNew), ranLevels = Gradient$rLNew, 
                                  studyDesign = Gradient$studyDesignNew, expected = expected_val, predictEtaMean = TRUE) |> 
              simplify2array()
            
            n_sites_batch   <- dim(predY_PA)[1]
            n_spp      <- dim(predY_PA)[2]
            n_post <- dim(predY_PA)[3]
            
            predY_hurdle <- array(
              0,
              dim = c(n_sites_batch, n_spp, n_post),
              dimnames = dimnames(predY_PA)
            )
            
            for (s in seq_len(n_post)) {
              
              th_s <- alfa_matrix[, s]
              is_present <- sweep(predY_PA[, , s], 2, th_s, ">=")
              
              # Apply threshold: keep probabilities where TRUE, zero out where FALSE, alfa
              pa_thresh <- predY_PA[, , s] * is_present
              
              # Construct hurdle: multiply thresholded probability by conditional abundance
              predY_hurdle[, , s] <- pa_thresh * predY_aCP[, , s]         
              
            }
            
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
            colnames(predY) <- colnames(hM_PA$Y)
              
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
              )
              
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
            # if(!dir.exists(part_path)) {
            # dir.create(part_path, recursive = TRUE, showWarnings = FALSE)
            # }
              
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
              writeLines(err_msg, file.path(pred_dir, paste0("ERROR_group_", sufix, group, "_batch_", b_num, ".txt")))
            })
          
          NULL
        }

stopCluster(cl)

message("    End time: ", Sys.time())