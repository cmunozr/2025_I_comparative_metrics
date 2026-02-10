#' DESCRIPTION:
#' This script calculates site-level biodiversity metrics (Taxonomic & Functional)
#' by comparing posterior predictions between matched Metso ( "Treatment") and 
#' Business as usual -BAU- ("Control") stands.
#' 
#' METHODOLOGY:
#' 1. **Batching:** Splits the m matches into manageable chunks (default: 500).
#' 2. **Lazy Loading:** Fetches posterior matrices from Arrow only for the active batch.
#' 3. **Parallelism:** Distributes the calculation of the active batch across N cores.
#' 4. **Storage:** Writes results to disk (Parquet) immediately after each batch
#' 
#' INPUTS:
#' - HMSC Model Metadata: For species names (`unfitted_...RData`).
#' - Traits Data: Functional traits for species (`TrData_complete.rds`).
#' - Predictions: Arrow dataset containing posterior predictions (`results/predictions`).
#' - Matched Pairs: Table of paired sites to compare (`matched_pairs_utm.rds`).
#' 
#' OUTPUTS:
#' - Directory: `results/metrics/[expected_type]/`
#' - Files: `metrics_batch_[N].parquet` (Long-format dataframe with metrics per posterior).
#' 
#' SYSTEM REQUIREMENTS:
#' - OS: Linux/macOS recommended (for efficient forking via `mclapply`).
#' - RAM: Depends on batch size (Batch 500 ~ requires approx 16-32GB RAM).

library(tidyverse)
library(Hmsc)
library(parallel)
library(here)
library(arrow)
library(FD)

os <- Sys.info()['sysname']

# If you are on a server (Linux), use 20. If Windows, usually stuck with 1 for mclapply
# or requires makeCluster for more. Assuming Linux/Server for 20 cores:
n_cores <- if (os == "Windows") 1 else 20

# Batch size: Control RAM usage. 
# 1000 rows * (posterior samples + species columns) = manageable chunk
batch_size <- 500

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

source(file.path("code", "_utilities_diversity_metrics_functions.R"))
source(file.path("code", "config_model.R"))
name_model <- generate_run_name(run_config)

expected_string <- "_expected_true" # or logic based on config

# 1. Data loading

# Load model metadata (Species names)
load(file.path("models", paste0("unfitted_", run_config$model_id, ".RData")))
name_spp <- colnames(models$fbs_M008$Y)

output_dir <- file.path("results", "metrics", expected_string)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

set.seed(11072024)

# Load Trait Data (Run once logic)
tr_proc_file <- file.path("results", "metrics", "TrData_processed.rds")

if (!file.exists(tr_proc_file)) {
  cat("----- Calculating TrData processing...\n")
  TrData <- readRDS(file.path("data", "traits", "TrData_complete.rds")) 
  rownames(TrData) <- name_spp
  TrData_numeric <- TrData |> 
    select(where(is.numeric), -c("Total.individuals", "Female", "Male", "Unknown", "Complete.measures"))
  
  dummy_a <- matrix(1, nrow = 1, ncol = length(name_spp))
  colnames(dummy_a) <- name_spp
  
  TrData_processed <- dbFD_preprocess_traits(x = TrData_numeric, a = dummy_a)
  saveRDS(TrData_processed, tr_proc_file)
} else {
  TrData_processed <- readRDS(tr_proc_file)
}

# Reduce PCoA space
eigs <- TrData_processed$eigenvalues
cum_var <- cumsum(eigs) / sum(eigs) * 100
index_95 <- cum_var <= 95
TrData_processed$traits.FRic <- TrData_processed$traits.FRic[ , index_95]

# 2. Prepare match and arrow connection

# Lazy connection to Arrow (Does not load data yet)
predY_ds <- open_dataset(
  file.path("results", "predictions"), 
  hive_style = TRUE, 
  unify_schemas = TRUE) |> 
  filter(
    modelid == run_config$model_id, 
    expected_type == expected_string, 
    year == 2021) 

# Get Control IDs
pred_id_control <- predY_ds |> 
  select(scenario, standid) |> 
  filter(scenario == "control") |> 
  collect() |> 
  pull(standid) |> 
  unique()

# Get Matches
matches <- readRDS(file.path("data", "metso", "matched_pairs_utm.rds")) |> 
  mutate(in_pred_control = standid_matched_control %in% pred_id_control, row_id = row_number()) |> 
  filter(in_pred_control == 1) 

matches <- matches[1:1000, ]

cat("Total matches to process:", nrow(matches), "\n")

# 3. Worker function to process a match pair

process_single_match <- function(idx, matches_df, pred_data_chunk, alpha, traits) {
  
  # 1. Get IDs for this specific match
  row_data <- matches_df[idx, ]
  ids_to_get <- c(row_data$standid_treated, row_data$standid_matched_control)
  
  # 2. Filter the pre-loaded chunk
  # Using base R subsetting for speed in parallel workers
  pred_sub <- pred_data_chunk[pred_data_chunk$standid %in% ids_to_get, ]
  
  # 3. Separate Matrices
  # Metso
  p_metso <- pred_sub[pred_sub$scenario == "metso", ]
  p_metso <- p_metso[order(p_metso$posterior), ] # Ensure order
  
  # Control == bau
  p_bau <- pred_sub[pred_sub$scenario == "control", ]
  p_bau <- p_bau[order(p_bau$posterior), ]
  
  # 4. Convert to Matrix (Drop non-species columns)
  # Adjust columns to drop based on your actual data structure
  cols_drop <- c("scenario", "standid", "posterior") 
  
  # Safety check, 0 posteriors in match
  if(nrow(p_metso) == 0 || nrow(p_bau) == 0) return(NULL)
  
  mat_metso <- as.matrix(p_metso[ , !(names(p_metso) %in% cols_drop)])
  mat_bau   <- as.matrix(p_bau[ , !(names(p_bau) %in% cols_drop)])
  
  # 5. Run Metrics
  if(nrow(mat_metso) == nrow(mat_bau)) {
    res_df <- tryCatch({
      # this function works only with expected = T predictions, in case of need expected = F it is needed another version
      # rest of the code can be used as it is
      calculate_metrics_vectorized(
        predY_metso_mat = mat_metso, 
        predY_bau_mat = mat_bau, 
        alpha = alpha,
        Traits = traits
      )
    }, error = function(e) return(NULL))
  } else {
    return(NULL)
  }
  
  if(!is.null(res_df)) {
    res_df$standid <- row_data$standid_treated
    res_df$match_id <- row_data$row_id
    res_df <- res_df[, c("standid", "match_id", "posterior", "metrics", "val")]
  }
}

# 4. Batch processing loop

# Split
indices <- seq_len(nrow(matches))
batches <- split(indices, ceiling(seq_along(indices)/batch_size))

start_time <- Sys.time()

cat(sprintf("[%s] Starting processing: %d batches | Cores: %d \n", 
            Sys.time(), length(batches), n_cores))

for(b_id in seq_along(batches)) {
  
  current_idx <- batches[[b_id]]
  cat(sprintf("[%s] --- Batch %d/%d (Rows %d-%d) ---\n", 
              Sys.time(), b_id, length(batches), min(current_idx), max(current_idx)))
  
  # A. data for this batch
  batch_matches <- matches[current_idx, ]
  needed_stands <- unique(c(batch_matches$standid_treated, batch_matches$standid_matched_control))
  needed_utm    <- unique(c(batch_matches$UTM200_metso, batch_matches$UTM200_control))
  
  # B. Collect data from Arrow 
  predY_chunk <- predY_ds |> 
    filter(utm_zone %in% needed_utm) |> 
    filter(standid %in% needed_stands) |> 
    select(scenario, standid, posterior, all_of(name_spp)) |> 
    collect()
  
  # C. Parallel Execution
  # mclapply uses forking on Linux, very memory efficient
  results <- mclapply(
    X = 1:nrow(batch_matches),
    FUN = function(i) {
      process_single_match(
        idx = i, 
        matches_df = batch_matches, 
        pred_data_chunk = predY_chunk, 
        alpha = 0.05, 
        traits = TrData_processed
      )
    },
    mc.cores = n_cores
  )
  
  # D. Save Batch
  results <- data.table::rbindlist(results)
  write_parquet(results, file.path(output_dir, paste0("metrics_batch_", b_id, ".parquet")))
  
  # E. Clean Memory
  rm(predY_chunk, results)
  gc() # Force garbage collection
}

end_time <- Sys.time()
total_time <- end_time - start_time

cat(sprintf("[%s] Done. Total time: %s\n", 
            Sys.time(), 
            format(total_time, digits = 2)))