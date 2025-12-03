# temporal for example to BES

library(sf)
library(dplyr)
library(here)

source(file.path("code", "config_model.R"))

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# Turn off S2 geometry to avoid spherical geometry errors with EPSG:4326
sf::sf_use_s2(FALSE)
set.seed(11072024)

# 1. SETUP & DATA LOADING
cat(sprintf("[%s] Starting setup...\n", Sys.time()))
sufix <- "control"

# Load sampling metso in S08a

samp <- readRDS(file.path("results", "predictions", "metso", "samp.rds"))

# Load Covariates
data_path <- file.path("data", "covariates", paste0("XData_hmsc_metso_", run_config$model_id, ".rds"))
XData_list <- readRDS(data_path)

samp <- sample(seq_len(length(XData_list$polygon_id)), size = 2000)

# Load Spatial Data
sp_df <- read_sf(file.path("data", "metso", "treatment_control_stand_filtered.gpkg"))

# Metso, replicating S08a
sp_df_metso <- sp_df |> 
  dplyr::filter(metso == 1, standid %in% XData_list$polygon_id) |>
  dplyr::distinct(standid, .keep_all = TRUE) |>
  dplyr::arrange(factor(standid, levels = XData_list$polygon_id))
sp_df_metso <- sp_df_metso[samp,]

# Control
sp_df_control <- sp_df |> 
  dplyr::filter(metso == 0) |>
  dplyr::distinct(standid, .keep_all = TRUE)

metso_samp <- sp_df_metso
matches <- data.frame("metso_standid" = NA, "control_standid" = NA)
matches.sp <- sp_df_control[0, ]

for(i in 1:nrow(metso_samp)){
  #i <- 17
  message(paste(i, "of", nrow(metso_samp)))
  metso_samp_i <- metso_samp[i, ]
  matches[i, "metso_standid"] <- metso_samp_i$standid
  
  metso_samp_buf_1km <- metso_samp_i |> 
    st_buffer(dist = 1000, endCapStyle = "ROUND") |> 
    summarise()
  metso_samp_buf_10km <- metso_samp_i |> 
    st_buffer(dist = 10000, endCapStyle = "ROUND") |> 
    summarise()
  donut <- st_difference(metso_samp_buf_10km, metso_samp_buf_1km)
  
  candidates <- sp_df_control[donut, c("standid", "area")]
    
  area_target_20per <- metso_samp_i$area*0.20
  area_target <- c((metso_samp_i$area-area_target_20per), (metso_samp_i$area+area_target_20per))
    
  candidates_area <- candidates |> 
    filter(area >= area_target[1], area <= area_target[2])
  

  used_controls <- na.omit(matches$control_standid)
  available_candidates <- candidates_area |>
    filter(!standid %in% used_controls)
  if(nrow(available_candidates) > 0){
    selection <- slice_sample(available_candidates, n = 1)
    matches.sp <- bind_rows(matches.sp, selection)
    matches[i, "control_standid"] <- selection$standid
  } else {
    matches[i, "control_standid"] <- NA 
  }
}

saveRDS(matches, file.path("data", "metso", "donut_matches.rds"))
sf::write_sf(matches.sp, file.path("data", "metso", "donut_matches_sp.gpkg"), 
             delete_layer = T)
  


