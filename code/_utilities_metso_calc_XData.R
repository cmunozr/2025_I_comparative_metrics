library(sf)
library(openxlsx)
library(tidyverse)
library(here)

source(file.path("code", "config_model.R"))
source(file.path("code", "_utilities_transform_covariates.R"))

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# Define
folder_name <- file.path("data", "covariates") 
run_calculate_XData <- T
run_new_var_XData <- F
dict_covar <- read.csv(file.path("data", "covariates", "dictionary_covariates.csv"), sep = ";")
mapping_functions <- read.xlsx(file.path("data", "covariates", "mapping_covariate_functions.xlsx"), sheet = 1)
sp_df <- sf::read_sf(file.path("data", "metso", "treatment_control_stand_filtered.gpkg")) 
sufix <- "control"

if(sufix == "metso"){
  sp_df <- sp_df |> dplyr::filter(metso == 1)  
}else{
  matches <- readRDS(file.path("data", "metso", "donut_matches.rds")) |> 
    na.omit()
  sp_df <- sp_df |> 
    dplyr::filter(metso == 0, standid %in% matches$control_standid) 
}

path_rds <- file.path(folder_name, paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))

# Call the function
construct_hmsc_XData(
  folder_name = folder_name,
  run_calc = run_calculate_XData, 
  run_new = run_new_var_XData,
  dict_covar = dict_covar,
  mapping_funcs = mapping_functions, 
  ref_data = sp_df, 
  data_sufix = sufix
)

run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)

XData_metso <- readRDS(path_rds)

XData_metso$XData <- XData_metso$XData |> 
  dplyr::select(dplyr::all_of(names(hM$XData)))

saveRDS(XData_metso, file = path_rds)

proc.time()
