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
sufix <- "metso"

# Call the function
XData_metso <- update_hmsc_data(
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

XData_metso <- XData_metso |> 
  dplyr::select(dplyr::all_of(names(hM$XData)))

path <- file.path(folder_name, paste0("XData_hmsc_", sufix, ".rds"))
saveRDS(XData_metso, file = path)
