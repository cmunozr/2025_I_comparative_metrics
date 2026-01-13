library(terra)
library(sf)
library(exactextractr)
library(purrr)
library(stringr)
library(maps)
library(ggplot2)
library(tidyr)
library(here)
library(dplyr)
library(zip)

os <- Sys.info()['sysname']

if (os == "Windows") {
  output_base_path <- "E:"
} else if (os %in% c("Linux", "Darwin")) {
  output_base_path <- "/home/avesta/munozcs/Documents/data/covariates/raw/"
  #output_base_path <- "/media/pitm/My Passport/"
}

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

source(file.path("code","_utilities_transform_covariates.R"))

log_info <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

log_info(paste("Running on:", os))
log_info(paste("Output path set to:", output_base_path))

#---------------

dict_covar <- read.csv(file.path("data", "covariates", "dictionary_covariates.csv"), sep = ";") # External 

MS_NFI_years <- c(2009, 2011, 2013, 2015, 2017, 2019, 2021)

set.seed(11072024)

#--------------
# Prepare paths and organize information for each set of covariates

## Three high Europe. Global no tiles. https://glad.umd.edu/users/Potapov/Europe_TCH/Tree_Height/. See _utilities_download_covariates.R code

tree_high_europe <- find_root_folder(folder_name = "tree_high_Europe")
full_path <- list.files(tree_high_europe, "*.tif$", full.names = T, recursive = F)

tree_high_europe <-
  tibble(path = full_path) |>
  tidyr::extract(
    col = path,
    into = c("var", "year"),
    regex = ".*/(.+?)_(\\d{4}).*\\.tif$", # divide a string in two and arrange in two different columns
    convert = TRUE,
    remove = FALSE
  ) |>
  dplyr::mutate(
    dataset = "tree_high_eu",
    UTM200 = NA,
    UTM10 = NA
  ) |>
  dplyr::select(dataset, var, path, year, UTM200, UTM10) |> 
  arrange(year) |> 
  filter(year %in% MS_NFI_years)


## Luke (Multi Source National Forest Inventory) before was organized as UTM-200 tiles in https://kartta.luke.fi/opendata/valinta.html 
## Now is a single layer for year/thematic  in http://www.nic.funet.fi/index/geodata/luke/vmi/. See _utilities_download_covariates.R script

luke <- find_root_folder("luke")
full_path <- list.files(luke, "*.tif$", full.names = T, recursive = T)

luke <- tibble(path = full_path) |>
  mutate(
    dataset = "luke",
    var = str_remove(basename(path), ".tif"),
    year = basename(dirname(path)),
    UTM200 = NA,
    UTM10 = NA
  ) |> 
  dplyr::select(dataset, var, path, year, UTM200, UTM10)|> 
  mutate(year = as.numeric(year)) |> 
  arrange(year)|> 
  filter(year %in% MS_NFI_years)

## Latvus (Canopy model) UTM 10 tiles. https://avoin.metsakeskus.fi/aineistot/Latvusmalli/Karttalehti/. See _utilities_download_covariates.R code

latvus <- find_root_folder("latvus")
full_path <- list.files(latvus, "*.tif$", full.names = T, recursive = T)

latvus <- tibble(path = full_path) |> 
  mutate(
    dataset = "latvus",
    year = basename(dirname(path)),
    var = paste0(dataset, "_", year),
    UTM200 = NA,
    UTM10 = str_extract(basename(path), "(?<=_)[^_.]+(?=\\_)")
  ) |> 
  dplyr::select(dataset, var, path, year, UTM200, UTM10) |> 
  mutate(year = as.numeric(year)) |> 
  arrange(year)|> 
  filter(year %in% MS_NFI_years)


## Elevation data, 2019 10x10 m Now is a single layer for year/thematic  in https://www.nic.funet.fi/index/geodata/mml/dem10m/2019/. See _utilities_download_covariates.R script

dem <- find_root_folder("DEM_Finland")
full_path <- list.files(dem, "*.tif$", full.names = T, recursive = T)

dem <- expand_grid(
    year = MS_NFI_years,  
    path = full_path
  ) |> 
  mutate(
    dataset = "DEM",
    var = "Elevation",
    year = year,
    UTM200 = str_remove(basename(path), ".tif"),
    UTM10 = NA
  ) |> 
  dplyr::select(dataset, var, path, year, UTM200, UTM10)

## Climatic data, Nordic gridded temperature and precipitation data from 1961 to present derived from in-situ observations
## 1 x 1 km, https://cds.climate.copernicus.eu/datasets/insitu-gridded-observations-nordic?tab=overview

clim <- find_root_folder("nordic_climate")
full_path <- list.files(clim, "*.tif$", full.names = T, recursive = F)

clim <-
  tibble(path = full_path) |>
  tidyr::extract(
    col = path,
    into = c("var", "year"),
    regex = ".*/(.+?)_(\\d{4}).*\\.tif$", # divide a string in two and arrange in two different columns
    convert = TRUE,
    remove = FALSE
  ) |>
  dplyr::mutate(
    dataset = "clim",
    UTM200 = NA,
    UTM10 = NA
  ) |>
  dplyr::select(dataset, var, path, year, UTM200, UTM10) |> 
  arrange(year)


#--------------------

# Prepare sampling sites (buffers level 2) and METSO and NON_METSO stands and routes polygons

coords <- list.files(file.path("data", "fbs"), pattern = "route_sections_L2", full.names = T) |> # external
  lapply(st_read)

metso <- st_read(file.path("data", "metso", "treatment_control_stand_v2.gpkg")) # external

utm_10 <- st_read(file.path("data", "utm35_zones", "TM35_karttalehtijako.gpkg"), layer = "utm10") # external
utm_200 <- st_read(file.path("data", "utm35_zones", "TM35_karttalehtijako.gpkg"), layer = "utm200") # external

metso_utm_join <- metso |> 
  dplyr::select(standid, metso) |> 
  mutate(set = "metso") |> 
  st_transform(st_crs(utm_200)) |> 
  st_join(y = utm_200) |> 
  rename(UTM200 = lehtitunnus) |> 
  st_join(y = utm_10) |> 
  rename(UTM10 = lehtitunnus) |> 
  dplyr::select(standid, metso, set, UTM200, UTM10)

metso_utm_join$poly_id <- 1:nrow(metso_utm_join) # id

coords_utm_join <- lapply(X = coords, FUN = function(X){
    X.i <- X |> 
      dplyr::select(sampleUnit, vakio, year) |> 
      mutate(set = "coords") |> 
      st_transform(st_crs(utm_200)) |> 
      st_join(y = utm_200) |> 
      rename(UTM200 = lehtitunnus) |> 
      st_join(y = utm_10) |> 
      rename(UTM10 = lehtitunnus) |> 
      dplyr::select(sampleUnit, vakio, year, set, UTM200, UTM10) |> 
      st_buffer(dist = 150) 
    X.i$poly_id <- 1:nrow(X.i) #id
    return(X.i)
    }
  )

#--------------------

# Checking extent of latvus, probably it would be better before download !!!! SOMETHING TO LEARN HERE (MASSIVE DATASET 300 GB at least)
if(os == "Windows"){
  zones <- c(10)
  datasets <- c("latvus")

  for(i in 1:length(datasets)){
    
    utm_zone_number <- zones[i]
    dataset_name <- datasets[i]
    
    utm_col_name <- paste0("UTM", utm_zone_number)
    
    utm_zone_dataset <- get(paste0("utm_", utm_zone_number)) |>
      rename({{ utm_col_name }} := "lehtitunnus") |>
      full_join(get(dataset_name), by = utm_col_name)
    
    utm_zone_dataset_coords <- lapply(X = coords_utm_join, FUN = function(X){
      temp <- utm_zone_dataset |> 
        left_join(st_drop_geometry(dplyr::select(X, -year)), by = utm_col_name, relationship = "many-to-many") |>
        filter(poly_id > 0)
      return(temp)
    } 
    )
    
    utm_zone_dataset_coords <- bind_rows(utm_zone_dataset_coords) |> 
      mutate(year = as.factor(year))
    
    fin <- st_as_sf(maps::map(database = "world", regions = "finland", plot = FALSE, fill = TRUE))
    
    p <- utm_zone_dataset_coords |> 
      filter(!is.na(year)) |>
      ggplot() +
      geom_sf(data = fin) +
      geom_sf(aes(color = year)) +
      facet_wrap(vars(year), ncol = 4)
    
    print(p)
  }
}

# Conclusion: Latvus must be discarded as a time-variant covariate for the Hmsc model (cannot handle NA values in XData), coverage is uneven and scattered

#-----------------------

master_covariates <- bind_rows(tree_high_europe, luke, clim) #dem
# master_covariates <- dem
write.csv(master_covariates, file.path("data", "master_covariates_temp.csv"), row.names = F)
coords_utm_join <- bind_rows(coords_utm_join)

#-----------------------

# Actually Get covariates information

master_covariates <- master_covariates |>
  prepare_raster_paths() |>
  group_by(dataset) |>
  mutate(raster_crs = get_raster_crs(first(proc_path))) |>
  ungroup() |>  
  filter(!is.na(raster_crs))

list_of_groups <- master_covariates |> 
  group_by(dataset, var, year) |> 
  group_split()

run <- T

log_info("Starting raw value extraction...")
if(run){
  log_info("Extracting coords data...")
  coords_raw_data <- purrr::map(.x = list_of_groups, 
                                ~extract_raw_values(polygons_sf = coords_utm_join, 
                                                    raster_info_df = .x, 
                                                    root_output_dir = output_base_path, 
                                                    id_column = "sampleUnit",
                                                    only_paths = F))
  log_info("Extracting METSO data...")
  metso_raw_data <- purrr::map(list_of_groups, 
                               ~extract_raw_values(metso_utm_join, .x,
                                                   root_output_dir = output_base_path, 
                                                   id_column = "standid",
                                                   only_paths = F, year_metso = 2021))
}else{
  coords_raw_data <- purrr::map(list_of_groups, 
                                ~extract_raw_values(coords_utm_join, .x, 
                                                    root_output_dir = output_base_path, 
                                                    id_column = "sampleUnit", 
                                                    only_paths = T))
  metso_raw_data <- purrr::map(list_of_groups, 
                               ~extract_raw_values(metso_utm_join, .x, 
                                                   root_output_dir = output_base_path, 
                                                   id_column = "standid", 
                                                   only_paths = T))
}
log_info("Extraction complete.")

#--------------------

# Aggregate Tile Covariate Data for All Polygon Types

log_info("Aggregating Coords Covariates...")
aggregate_covariates(
  raw_data_list = coords_raw_data,
  polygon_type = coords_utm_join$set[1],
  dict_covar = dict_covar,
  output_root_dir = output_base_path,
  df = st_drop_geometry(coords_utm_join)
)

log_info("Aggregating METSO Covariates...")
aggregate_covariates(
  raw_data_list = metso_raw_data,
  polygon_type = metso_utm_join$set[1],
  dict_covar = dict_covar,
  output_root_dir = output_base_path,
  df = st_drop_geometry(metso_utm_join)
)

#-------------------

toKeep <- c("dict_covar", "output_base_path", "coords_utm_join", "metso_utm_join",
            "prepare_covariates_data_toplot", "generate_covariate_plot", "generate_binary_plot",
            "log_info", "os")
rm(list = setdiff(ls(), toKeep));gc()

# transform to year separated data to one consolidated for covariate for all polygon type

folder_name <- file.path("data", "covariates", "pre_processed") 

dir.create(folder_name, showWarnings = F)
  
covar_nm <- dict_covar |> 
  filter(processed == 1) |> 
  pull(var_shultz) |> 
  unique()

rds_data_dir <- file.path(output_base_path,"intermediate_aggregated")
  
polygon_types <- c(coords_utm_join$set[1], metso_utm_join$set[1])
#polygon_types <- c(coords_utm_join$set[1]) #just coords of biotopes

for(p in 1:length(polygon_types)){

  rds_files <- list.files(rds_data_dir, pattern = paste0(polygon_types[p], ".rds$"), recursive = TRUE, full.names = TRUE)
  
  for (i in 1:length(covar_nm)) {
    # i <- 1
    covar_nm_i <- covar_nm[i] 
    log_info(paste("Merging and saving:", i, "of", length(covar_nm), "variables -", covar_nm_i))
    
    relevant_files <- rds_files[str_detect(basename(rds_files), covar_nm_i)]
    
    if(length(relevant_files) == 0){
      next
    } 
    
    covar_data_i <- lapply(relevant_files, readRDS) |> 
      bind_rows()
    
    file_path <- file.path(folder_name,  paste0(covar_nm_i, "_", polygon_types[p], ".rds"))
    saveRDS(covar_data_i, file = file_path)
  }
}

#--------------------

toKeep <- c("covar_nm", "coords_utm_join", "metso_utm_join", "dict_covar", "folder_name",
            "prepare_covariates_data_toplot", "generate_covariate_plot" , "generate_binary_plot",
            "log_info", "os")
rm(list = setdiff(ls(), toKeep));gc()

if(os != "Windows"){
  stop("finishing routine for linux server")
}

# processing to plot sample

rds_files <- list.files(folder_name, pattern = ".rds$", recursive = TRUE, full.names = TRUE)
  
for (i in 1:length(covar_nm)) {
  
  covar_nm_i <- covar_nm[i] 
  
  message(paste("processing for plot", i, "of", length(covar_nm), "variables", "(", covar_nm_i, ")" ))
  
  relevant_files <- rds_files[str_detect(basename(rds_files), covar_nm_i)]
  
  if(length(relevant_files) == 0){
    next
  } 
  
  coords <- readRDS(str_subset(relevant_files, coords_utm_join$set[1]))
  
  n_ <- coords$polygon_id |> unique() |> length()
  
  metso_ <- readRDS(str_subset(relevant_files, metso_utm_join$set[1])) |> 
    mutate(polygon_id = as.integer(polygon_id))
  
  n_of_metso_stands_  <-  length(unique(metso_utm_join$standid))
  
  if(n_ > n_of_metso_stands_){
    n_ <- n_of_metso_stands_
  }
  
  ids_ <- metso_utm_join |>
    st_drop_geometry() |> # Drop the geometry column first
    mutate(poly_id = as.integer(poly_id)) |> 
    group_by(metso) |>
    slice_sample(n = n_, replace = FALSE) |>
    ungroup()
  
  metso_ <- left_join(metso_, ids_, join_by("polygon_id" == "poly_id")) |> 
    filter(!is.na(metso)) |> 
    mutate(polygon_id = as.factor(polygon_id))

  covar_data_i <- bind_rows(coords, metso_) |> 
    dplyr::select(colnames(coords), metso)
  
  rm(coords, metso_); gc()
  
  covar_data_i <- prepare_covariates_data_toplot(covar_data = covar_data_i)
  
  is_binary <- dict_covar[which(dict_covar$var_shultz == covar_data_i$variable[1]), "binary"] |> 
    as.logical()
  
  if (is_binary) {
    generate_binary_plot(covariate_name = covar_data_i$variable[1], data = covar_data_i, folder = folder_name)
  } else {
    generate_covariate_plot(covariate_name = covar_data_i$variable[1], data = covar_data_i, folder = folder_name)
  }  
}

