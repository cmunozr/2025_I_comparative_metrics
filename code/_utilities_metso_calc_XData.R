library(sf)
library(openxlsx)
library(dplyr)
library(here)
library(stringr)

source(file.path("code", "config_model.R"))
source(file.path("code", "_utilities_transform_covariates.R"))

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# 2. Define
folder_name <- file.path("data", "covariates") 
run_calculate_XData <- T
run_new_var_XData <- F
dict_covar <- read.csv(file.path("data", "covariates", "dictionary_covariates.csv"), sep = ";")
mapping_functions <- read.xlsx(file.path("data", "covariates", "mapping_covariate_functions.xlsx"), sheet = 1)
sp_df <- sf::read_sf(file.path("data", "metso", "treatment_control_stand_v4.gpkg")) 
trt.control <- read.csv(file.path("data", "metso", "treatment_control_standid.csv")) |> 
  na.omit()

# 3. to run the function construct_hmsc_XData the pre-processed variables need to be split between control and treated (metso)
# one time process

path_preprocessed <- list.files(file.path("data", "covariates", "pre_processed"), 
                                pattern = paste0("metso", ".rds$"), full.names = T, 
                                recursive = T)
pre_process <- F

if(pre_process){
  lapply(X = path_preprocessed, FUN = function(X){
    # X <- path_preprocessed[1]
    var <- readRDS(X) 
    c <- var |> 
      filter(polygon_id %in% trt.control$standid[trt.control$match_role == "control"])
    new_path <- sub(pattern = "metso", replacement = "control", x = X)
    saveRDS(object = c, file = new_path)
    m <- var |> 
      filter(polygon_id %in% trt.control$standid[trt.control$match_role == "treated"])
    saveRDS(object = m, file = X)
    invisible("success")
  })
}

# 4. each type of polygon needs to be run by itself

sufix <- "metso"

if(sufix == "metso"){

  sp_df <- sp_df |> dplyr::filter(metso == 1)  

}else if(sufix == "control"){
  
  sp_df <- sp_df |> dplyr::filter(metso == 0) 
}

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

# 5. control/treated(metso) covariates need to be the same as the fitted model   

write_complete_covariates <- F

run_name <- generate_run_name(run_config)
fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
hM <- readRDS(fitted_full_model_path)

if(write_complete_covariates){
  path_rds <- file.path(folder_name, paste0("XData_hmsc_", sufix, "_", run_config$model_id, "_complete.rds"))
}else{
  path_rds <- file.path(folder_name, paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
}

XData_ <- readRDS(path_rds)

XData_$XData <- XData_$XData |> 
  dplyr::select(dplyr::any_of(names(hM$XData)))

if(write_complete_covariates){
  path_rds <- file.path(folder_name, paste0("XData_hmsc_", sufix, "_", run_config$model_id, ".rds"))
}

saveRDS(XData_, file = path_rds)


