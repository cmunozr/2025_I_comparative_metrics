library(Hmsc)
library(sf)
library(here)

source(file.path("code", "config_model.R"))

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# 1. Define & Setup
cat("Starting setup...\n")
sufix <- "metso"

# Load Covariates
data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, ".rds"))
XData_metso <- readRDS(data_path) |> 
  as.data.frame()

# Load Spatial Data and Transform to METERS (EUREF-FIN)
# EPSG:3067 is crucial for Finland to keep units in meters for HMSC
sp_df <- sf::read_sf(file.path("data", "metso", "treatment_control_stand_filtered.gpkg")) |> 
  st_centroid() |> 
  st_transform(crs = "EPSG:4326") 

# Load Model
run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)

# Extract Spatial Random Level Name automatically
spatial_level_name <- hM$rLNames[1] 
cat(sprintf("Spatial random level identified as: %s\n", spatial_level_name))


# 2. AUTOMATED BATCH PREDICTION LOOP


# Configuration for batching
total_rows <- nrow(XData_metso)
batch_size <- 2000  # Adjust based on RAM. 2000 is usually safe.
num_batches <- ceiling(total_rows / batch_size)

# Create output directory for predictions if it doesn't exist
pred_dir <- file.path("results", "predictions", sufix)
if(!dir.exists(pred_dir)) dir.create(pred_dir, recursive = TRUE, showWarnings = F)

cat(sprintf("Starting prediction loop for %d stands in %d batches...\n", total_rows, num_batches))

for (i in 1:num_batches) {
  
  # A. Define indices for this batch
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, total_rows)
  indices <- start_idx:end_idx
  
  cat(sprintf("Processing Batch %d/%d (Rows %d-%d)... ", i, num_batches, start_idx, end_idx))
  
  # B. Check if file already exists (Skip if done)
  batch_file <- file.path(pred_dir, paste0("predY_", run_config$model_id, "_batch_", i, ".rds"))
  if (file.exists(batch_file)) {
    cat("Skipping (Already exists)\n")
    next
  }
  
  t_start <- proc.time()
  
  # C. Prepare Data Slices
  XData_sub <- XData_metso[indices, ]
  
  # Extract coordinates for this batch
  coords_sub <- st_coordinates(sp_df)[indices, ]
  
  # Create the Named List dynamically
  sDataNew_sub <- list()
  sDataNew_sub[[spatial_level_name]] <- coords_sub
  
  # D. Prepare Gradient

  try({
    Gradient <- prepareGradient(hM, XDataNew = XData_sub, sDataNew = sDataNew_sub)
    predY <- predict(hM, Gradient = Gradient, expected = FALSE)
    saveRDS(predY, file = batch_file)
    rm(predY, Gradient, XData_sub, sDataNew_sub); gc()
  })
  
  t_end <- proc.time() - t_start
  cat(sprintf("Done in %.1f sec.\n", t_end[3]))
}
cat("-----------DONE---------------/n")
