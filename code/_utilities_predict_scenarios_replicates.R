# it replicates n times the prediction over a handful of standids, to check how vary the data in different circunstances
# of creation, it works with the n-dimensional array

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
sampling <- T

# 1. SETUP & DATA LOADING
cat(sprintf("[%s] Starting setup...\n", Sys.time()))
sufix <- "metso" # "metso" <------------- Change
expected_val <- T

if(expected_val){
  expected_string <- "_expected_true"
}else{
  expected_string <- "_expected_false"
}

# Create output directory

pred_dir <- file.path("results", paste0("predictions_", expected_string, "_replicates"))
pred_dir_sufix <- file.path(pred_dir, sufix)

if(!dir.exists(pred_dir_sufix)) dir.create(pred_dir_sufix, recursive = TRUE, showWarnings = FALSE)

# Call spatial data frame (this df has filter check the metso setup)
sp_df <- read_sf(here("data", "metso", "treatment_control_stand_v2.gpkg"))

# call matches 
matches <- readRDS(file.path("data", "metso", "raw", "matched_pairs.rds")) |> 
  na.omit()

# call XData

data_path <- file.path("data", "covariates", paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)
XData <- XData_list$XData |> 
  as.data.frame()

# in case of sampling for test u other reason

if(sampling){
  
  matches <- matches |> 
    filter(standid_treated %in% sp_df$standid | standid_matched_control %in% sp_df$standid)
  
  if(sufix == "metso"){
    samp <- sample(seq_len(nrow(matches)), size = 50)
    matches <- matches[samp,]
  }
  
  if(sufix == "control"){
    samp_id <- read_rds(file.path(pred_dir, "metso", "pred_ids.rds"))  
    matches <- matches |> 
      filter(standid_treated %in% samp_id)
  }
} 

if(sufix == "metso"){
  # Load Spatial Data
  
  sp_df <- sp_df |>
    st_centroid() |>
    st_transform("EPSG:4326") |>
    dplyr::filter(metso == 1, standid %in% matches$standid_treated) |>
    dplyr::distinct(standid, .keep_all = TRUE) |>
    dplyr::arrange(factor(standid, levels = XData_list$polygon_id))
  
}else if(sufix == "control"){
  sp_df <- sp_df |> 
    st_centroid() |>
    st_transform("EPSG:4326") |>
    dplyr::filter(metso == 0, standid %in% matches$standid_matched_control) |>
    dplyr::distinct(standid, .keep_all = TRUE) |>
    dplyr::arrange(factor(standid, levels = XData_list$polygon_id))
}

index <- XData_list$polygon_id %in% sp_df$standid
XData <- XData[index, ]
XData_list$polygon_id <- XData_list$polygon_id[index]
sp_df <- sp_df[sp_df$standid %in% XData_list$polygon_id, ]

saveRDS(sp_df$standid, file.path(pred_dir_sufix, "pred_ids.rds"))

coords <- st_coordinates(sp_df) |> 
  as.data.frame()


# Ensure alignment between XData and sp_df after filtering
match_idx <- match(sp_df$standid, XData_list$polygon_id)
XData <- XData[match_idx, , drop = FALSE]
stopifnot(all(XData_list$polygon_id[match_idx] == sp_df$standid))


# Load Model
run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)

# Extract Spatial Random Level Name
# This ensures the list passed to prepareGradient is named correctly (e.g., "route" or "Route")
spatial_level_name <- hM$rLNames[1] 
cat(sprintf("Spatial random level identified as: %s\n", spatial_level_name))

# C. Prepare Data Slices
# Slice covariates
XData_sub <- XData
        
# Slice coordinates
coords_sub <- coords
              
# Create Named List for Gradient (Crucial for prepareGradient)
sDataNew_sub <- list()
sDataNew_sub[[spatial_level_name]] <- coords_sub


# ---------------------------
# PARALLEL REPLICATES
# ---------------------------

library(doParallel)
library(foreach)

n_replicates <- 10
n_cores <- 10  # You can use 10 workers for 10 replicates

cat(sprintf("Initializing parallel cluster with %d workers...\n", n_cores))
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)

parallel_time <- system.time({
  
  foreach(i = 1:n_replicates,
          .packages = c("Hmsc", "sf"),
          .export = c("XData_sub", "coords_sub", "sDataNew_sub",
                      "hM", "run_config", "spatial_level_name",
                      "pred_dir_sufix", "expected_val")) %dopar% {
                        
                        cat(sprintf("[%s] Worker %d running replicate %d\n",
                                    Sys.time(), Sys.getpid(), i))
                        
                        batch_filename <- paste0("predY_iteration",
                                                 run_config$model_id,
                                                 "_replicate_", i, ".rds")
                        batch_file_path <- file.path(pred_dir_sufix, batch_filename)
                        
                        tryCatch({
                          
                          # Prepare gradient
                          Gradient <- prepareGradient(
                            hM,
                            XDataNew = XData_sub,
                            sDataNew = sDataNew_sub
                          )
                          
                          # Predict deterministically
                          predY <- predict(
                            hM,
                            Gradient = Gradient,
                            expected = expected_val,
                            predictEtaMean = TRUE
                          )
                          
                          saveRDS(predY, batch_file_path)
                          
                          rm(predY, Gradient)
                          gc()
                          
                          cat(sprintf("[%s] Worker %d DONE replicate %d\n",
                                      Sys.time(), Sys.getpid(), i))
                          
                        }, error = function(e) {
                          msg <- paste("Error replicate", i, ":", conditionMessage(e))
                          cat(msg, "\n")
                          writeLines(msg,
                                     file.path(pred_dir_sufix,
                                               paste0("ERROR_replicate_", i, ".txt")))
                        })
                        
                        NULL
                      }
})

stopCluster(cl)

cat(sprintf("Finished parallel replicates in %.2f minutes\n",
            parallel_time["elapsed"] / 60))

#---------------------------
# CHECKING

rm(list = ls()); gc()

prediction_files <- list.files(file.path("results", paste0("predictions", "_expected_true", "_replicates"), "metso"), pattern = "predY", full.names = TRUE)

predY <- lapply(X = prediction_files, FUN = function(X){readRDS(X) |> 
  simplify2array()})

dims_ok <- vapply(
  predY, 
  function(a) paste(dim(a), collapse = "x"), 
  character(1)
)

predY_4d <- simplify2array(predY)

dim(predY_4d)

# sites, species, posterior sample, replicate
predY_4d[,8,10,1] %>% summary()
predY_4d[,8,11,2] %>% summary()
predY_4d[,8,12,4] %>% summary()
predY_4d[,8,13,5] %>% summary()
predY_4d[,8,14,6] %>% summary()
predY_4d[,8,15,7] %>% summary()

   