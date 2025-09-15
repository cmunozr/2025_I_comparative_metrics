library(terra)
library(sf)
library(exactextractr)
library(purrr)
library(stringr)
library(maps)
library(ggplot2)
library(tidyr)
library(dplyr)

source("R/_utilities_transform_covariates.R")
dict_covar <- read.csv("data/covariates/dictionary_covariates.csv", sep = ";") # External 

set.seed(11072024)
select
#--------------
# Prepare paths and organize information for each set of covariates

## Three high Europe. Global no tiles. https://glad.umd.edu/users/Potapov/Europe_TCH/Tree_Height/. See _utilities_download_covariates.R code

tree_high_europe <- find_root_folder("tree_high_Europe")
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
  filter(year %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021))


## Luke (Multi Source National Forest Inventory) before UTM-200 tiles. http://www.nic.funet.fi/index/geodata/luke/vmi/. See _utilities_download_covariates.R code

luke <- find_root_folder("luke_")
full_path <- list.files(luke, "*.tif$", full.names = T, recursive = F)

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
  filter(year %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021))

## Latvus (Canopy model) UTM 10 tiles. https://avoin.metsakeskus.fi/aineistot/Latvusmalli/Karttalehti/. See _utilities_download_covariates.R code

latvus <- find_root_folder("latvus")
full_path <- list.files(latvus, "*.tif$", full.names = T, recursive = F)

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
  filter(year %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021))

#--------------------

# Prepare sampling sites (buffers level 2) and METSO and NON_METSO stands and routes polygons

coords <- list.files(file.path("data", "fbs"), pattern = "route_sections_L2", full.names = T) |> # external
  lapply(st_read)

metso <- st_read("data/metso/treatment_control_stand.gpkg") # external

utm_10 <- st_read("data/utm35_zones/TM35_karttalehtijako.gpkg", layer = "utm10") # external
utm_200 <- st_read("data/utm35_zones/TM35_karttalehtijako.gpkg", layer = "utm200") # external

metso_utm_join <- metso |> 
  dplyr::select(standid, metso) |> 
  mutate(set = "metso") |> 
  st_transform(st_crs(utm_200)) |> 
  st_join(y = utm_200) |> 
  rename(UTM200 = lehtitunnus) |> 
  st_join(y = utm_10) |> 
  rename(UTM10 = lehtitunnus) |> 
  dplyr::select(standid, metso, set, UTM200, UTM10)

sampling_metso <- T

if(sampling_metso){
  metso_utm_join <- metso_utm_join |> 
    group_by(metso) |>
    sample_n(size = 10000, replace = FALSE) |>
    ungroup()  
}

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

# Conclusion: Latvus must be discarded as a time-variant covariate for the Hmsc model (cannot handle NA values in XData), coverage is uneven and scattered

#-----------------------

master_covariates <- bind_rows(tree_high_europe, luke)
write.csv(master_covariates, "data/master_covariates_temp.csv", row.names = F)

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

if(run){
  coords_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(coords_utm_join, .x, id_column = "sampleUnit"))
  metso_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(metso_utm_join, .x, id_column = "standid"))
}else{
  coords_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(coords_utm_join, .x, id_column = "sampleUnit", only_paths = T))
  metso_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(metso_utm_join, .x, id_column = "standid", only_paths = T))
}

toKeep <- c("coords_raw_data", "metso_raw_data", "utm_10", "utm_200", "master_covariates")
rm(list = setdiff(ls(), toKeep));gc()

#--------------------

# Aggregate Tile Covariate Data for All Polygon Types

root_output_dir <- file.path("D:", "intermediate_aggregated")

aggregate_tile_covariates(
  raw_data_list = coords_raw_data,
  polygon_type = coords_utm_join$set[1],
  dict_covar = dict_covar,
  output_root_dir = root_output_dir,
  df = st_drop_geometry(coords_utm_join)
)

aggregate_covariates( # havent run yet
  raw_data_list = metso_raw_data,
  polygon_type = metso_utm_join$set[1],
  dict_covar = dict_covar,
  output_root_dir = root_output_dir,
  df = st_drop_geometry(metso_utm_join)
)