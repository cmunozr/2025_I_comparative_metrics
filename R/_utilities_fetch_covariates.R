library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

source("R/_utilities_transform_covariates.R")

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


## Luke (Multi Source National Forest Inventory) UTM-200 tiles. https://kartta.luke.fi/opendata/valinta.html. See _utilities_download_covariates.py code

luke <- find_root_folder("luke")
full_path <- list.files(luke, "*.zip$", full.names = T, recursive = F)

# Check and split if the utm code actually has multiple utm 200 zones, example: KL2 
split_utm_code <- function(code) {
  if (str_detect(code, "[A-Z]{2,}")) {
    letters_part <- str_extract(code, "[A-Za-z]+")
    number_part <- str_extract(code, "\\d+")
    individual_letters <- str_split_1(letters_part, "")
    return(paste0(individual_letters, number_part))
  } else {
    return(code)
  }
}

luke <- tibble(path = full_path) |>
  mutate(
    dataset = "luke",
    var = str_extract(basename(path), "^(.*)_[^_]*$", group = 1),
    year = basename(dirname(path)),
    UTM200 = str_extract(basename(path), "(?<=_)[^_.]+(?=\\.)"),
    UTM10 = NA
  ) |> 
  mutate(UTM200 = purrr::map(UTM200, split_utm_code)) |> 
  unnest(UTM200) |> 
  select(dataset, var, path, year, UTM200, UTM10)|> 
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
  select(dataset, var, path, year, UTM200, UTM10) |> 
  mutate(year = as.numeric(year)) |> 
  arrange(year)|> 
  filter(year %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021))

#--------------------

# Prepare sampling sites (buffers level 2) and METSO and NON_METSO stands and routes polygons

coords <- list.files(file.path("data", "fbs"), pattern = "route_sections_L2", full.names = T) |> 
  lapply(st_read)

metso <- st_read("data/metso/treatment_control_stand.gpkg")

utm_10 <- st_read("data/utm35_zones/TM35_karttalehtijako.gpkg", layer = "utm10")
utm_200 <- st_read("data/utm35_zones/TM35_karttalehtijako.gpkg", layer = "utm200")

metso_utm_join <- metso |> 
  select(standid, metso) |> 
  mutate(set = "metso") |> 
  st_transform(st_crs(utm_200)) |> 
  st_join(y = utm_200) |> 
  rename(UTM200 = lehtitunnus) |> 
  st_join(y = utm_10) |> 
  rename(UTM10 = lehtitunnus) |> 
  select(standid, metso, set, UTM200, UTM10)

metso_utm_join$poly_id <- 1:nrow(metso_utm_join) # id

coords_utm_join <- lapply(X = coords, FUN = function(X){
    X.i <- X |> 
      select(sampleUnit, vakio, year) |> 
      mutate(set = "coords") |> 
      st_transform(st_crs(utm_200)) |> 
      st_join(y = utm_200) |> 
      rename(UTM200 = lehtitunnus) |> 
      st_join(y = utm_10) |> 
      rename(UTM10 = lehtitunnus) |> 
      select(sampleUnit, vakio, year, set, UTM200, UTM10) |> 
      st_buffer(dist = 150) 
    X.i$poly_id <- 1:nrow(X.i) #id
    return(X.i)
    }
  )

#--------------------

# Checking extent of latvus, probably it would be better before download !!!! SOMETHING TO LEARN HERE (MASSIVE DATASET 300 GB at least)

zones <- c(10, 200)
datasets <- c("latvus", "luke")

for(i in 1:length(datasets)){
  
  utm_zone_number <- zones[i]
  dataset_name <- datasets[i]
  
  utm_col_name <- paste0("UTM", utm_zone_number)
  
  utm_zone_dataset <- get(paste0("utm_", utm_zone_number)) |>
    rename({{ utm_col_name }} := "lehtitunnus") |>
    full_join(get(dataset_name), by = utm_col_name)
  
  utm_zone_dataset_coords <- lapply(X = coords_utm_join, FUN = function(X){
    temp <- utm_zone_dataset |> 
      left_join(st_drop_geometry(select(X, -year)), by = utm_col_name, relationship = "many-to-many") |>
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

coords_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(coords_utm_join, .x, id_column = "sampleUnit"))
metso_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(metso_utm_join, .x, id_column = "standid"))

toKeep <- c("coords_raw_data", "metso_raw_data", "utm_10", "utm_200", "master_covariates")
rm(list = setdiff(ls(), toKeep));gc()

#--------------------

# Agregate to polygon level

## because problems in variable names, I needed to construct a dictionary

dict_covar <- read.csv("data/covariates/dictionary_covariates.csv", sep = ";")
root_name <- file.path("D:", "intermediate_aggregated")
dir.create(root_name, showWarnings = F)
polygon_type <- "metso" #  "metso" = metso_raw_data, "coords" = coords_raw_data

folder_name <- file.path(root_name, polygon_type)
dir.create(folder_name, showWarnings = F)

for(a in 1:length(metso_raw_data)){ #<---------- REPLACE
  #a <- 5
  message(paste("running", a, "of", length(metso_raw_data))) #<---------- REPLACE
  dfa <- purrr::map_dfr(.x = metso_raw_data[[a]], #<---------- REPLACE
                      .f = ~ compact(readRDS(str_replace(.x, pattern = "results", replacement = "D:")))
  )
  if(!is.na(dfa$variable[1])){
    index <- which(dfa$variable[1] == dict_covar$var_error)
    if(dfa$dataset[1] == "luke"){
      dfa$variable <- dict_covar$var_new[index]  
    }else{
      dfa$variable <- dfa$variable
    }
  }else{
    dfa$variable <- "Tree_removal_latest_date"
  }
  
  dfa$polygon_source <- polygon_type
  
  file_name <-  file.path(folder_name, paste0(dfa$variable[1], "_", dfa$year[1], "_", polygon_type, ".rds"))
  saveRDS(object = dfa, file = file_name)
  rm(dfa); gc()
}
