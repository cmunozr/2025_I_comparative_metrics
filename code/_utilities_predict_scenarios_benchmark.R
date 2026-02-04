cat("-----------START PARALLEL BENCHMARKING--------------\n")
library(Hmsc)
library(sf)
library(here)
library(doParallel)
library(foreach)
library(tidyverse)

if(!(Sys.getenv("RSTUDIO") == "1")) {
  setwd(here::here()) 
}

source(file.path("code", "config_model.R"))

# 1. CONFIGURATION
sufix <- "control" 
expected_val <- TRUE
test_sizes <- c(50, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000) 
n_cores <- length(test_sizes) # Uses 8 cores (1 per test case)

results_file <- file.path("results", "benchmark_parallel_results.csv")

# 2. SETUP & DATA LOADING

sf::sf_use_s2(FALSE)
set.seed(11072024)

cat("Loading Data...\n")
sp_df <- read_sf(here("data", "metso", "treatment_control_stand_v2.gpkg"))
matches <- readRDS(file.path("data", "metso", "raw", "matched_pairs.rds"))

data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)
XData <- as.data.frame(XData_list$XData)

# Filter & Align
target_ids <- matches$standid_matched_control
sp_df <- sp_df |>
  dplyr::filter(metso == 0, standid %in% target_ids) |>
  dplyr::distinct(standid, .keep_all = TRUE) |> 
  st_centroid() |>
  st_transform("EPSG:4326")

valid_ids <- intersect(sp_df$standid, XData_list$polygon_id)
sp_df <- sp_df[sp_df$standid %in% valid_ids, ]
match_idx <- match(sp_df$standid, XData_list$polygon_id)
XData <- XData[match_idx, , drop = FALSE]
coords <- st_coordinates(sp_df) |> as.data.frame()

# Subset Data (We load enough rows for the largest batch to run)
# Since we are not looping the data, we just need enough rows to satisfy the largest request (1000)
max_needed <- max(test_sizes)
index <- sample.int(nrow(coords), max_needed)
XData_bench <- XData[index, ]
coords_bench <- coords[index, ]

# Load Model
run_name <- generate_run_name(run_config)
hM <- readRDS(file.path("models", run_name, paste0("fitted_", run_name, ".rds")))
spatial_level_name <- hM$rLNames[1] 

cat(sprintf("Data ready. Launching %d parallel workers on %d cores.\n", length(test_sizes), n_cores))

# 3. INITIALIZE CSV
# Write headers immediately
headers <- data.frame(BatchSize=integer(), TotalTime_Sec=numeric(), Sec_Per_Row=numeric(), Status=character())
write.csv(headers, results_file, row.names = FALSE)
cat(sprintf("Writing live results to: %s\n", results_file))

# 4. PARALLEL EXECUTION
cl <- makeCluster(n_cores, outfile = "") 
registerDoParallel(cl)

# Foreach iterates through the vector of sizes
# Each size is sent to a different worker
out <- foreach(size = test_sizes, 
               .packages = c("Hmsc"), 
               .export = c("XData_bench", "coords_bench", "hM", "spatial_level_name", "expected_val", "results_file"),
               .errorhandling = "pass") %dopar% {
                 
                 cat(sprintf("Worker %d STARTING batch size: %d\n", Sys.getpid(), size))
                 
                 # We use tryCatch to ensure even if a batch fails, we record the error to CSV
                 status <- "SUCCESS"
                 timer <- system.time({
                   tryCatch({
                     # Slice the exact number of rows for this specific test
                     XData_sub <- XData_bench[1:size, ]
                     coords_sub <- coords_bench[1:size, ]
                     
                     sDataNew_sub <- list()
                     sDataNew_sub[[spatial_level_name]] <- coords_sub
                     
                     # Predict
                     Gradient <- prepareGradient(hM, XDataNew = XData_sub, sDataNew = sDataNew_sub)
                     predY <- predict(hM, Gradient = Gradient, expected = expected_val, predictEtaMean = TRUE)
                     
                   }, error = function(e) {
                     status <<- paste("ERROR:", e$message)
                   })
                 })
                 
                 # Prepare Result Row
                 elapsed <- timer["elapsed"]
                 sec_row <- elapsed / size
                 
                 # Create data frame 
                 res_row <- data.frame(
                   BatchSize = size,
                   TotalTime_Sec = elapsed,
                   Sec_Per_Row = sec_row,
                   Status = status
                 )
                 
                 # WRITE TO CSV IMMEDIATELY (Append mode)
                 write.table(res_row, file = results_file, append = TRUE, sep = ",", 
                             col.names = FALSE, row.names = FALSE, quote = TRUE)
                 
                 cat(sprintf("Worker %d FINISHED batch %d in %.2fs (%.4f s/row)\n", Sys.getpid(), size, elapsed, sec_row))
                 
                 # Return NULL to save memory on master
                 NULL
               }

stopCluster(cl)
cat("-----------DONE---------------\n")

#------------------
library(ggplot2)

rm(list = ls()); gc()

time_batches <- read.csv(file.path("results", "benchmark_parallel_results_.csv"), sep =";")
time_batches$etaMean <- c(rep(T, 7), F, rep(T, 5), rep(F, 4))

time_batches |> 
  ggplot(mapping = aes(x = BatchSize, y = Sec_Per_Row, colour = etaMean)) +
  geom_point() +
  xlim(c(0,500))