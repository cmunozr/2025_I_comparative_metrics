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
sp_df <- sf::read_sf(file.path("data", "metso", "treatment_control_stand_v2.gpkg")) 
matches <- readRDS(file.path("data", "metso", "raw", "matched_pairs.rds")) |> 
  na.omit()
sufix <- "control"

if(sufix == "metso"){

  sp_df <- sp_df |> dplyr::filter(metso == 1)  

}else if(sufix == "control"){
  
  path_preprocessed <- list.files(file.path("data", "covariates", "pre_processed"), 
                                  pattern = paste0("metso", ".rds$"), full.names = T, 
                                  recursive = T)
  lapply(X = path_preprocessed, FUN = function(X){
    var <- readRDS(X) 
    c <- var |> 
      filter(polygon_id %in% matches$standid_matched_control)
    new_path <- sub(pattern = "metso", replacement = "control", x = X)
    saveRDS(object = c, file = new_path)
    m <- var |> 
      filter(polygon_id %in% matches$standid_treated)
    saveRDS(object = m, file = X)
    invisible("success")
  })
  
  sp_df <- sp_df |> 
    dplyr::filter(metso == 0, standid %in% matches$standid_matched_control) 
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

XData_ <- readRDS(path_rds)

XData_$XData <- XData_$XData |> 
  dplyr::select(dplyr::all_of(names(hM$XData)))

saveRDS(XData_, file = path_rds)

proc.time()
