library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(purrr)

# Prepare paths and organize information for each set of covariates

## Get root folder to creathe the paths for a set of covariates, example: "latvus". 
## It iterates over all the hard disks available

find_root_folder <- function(folder_name) {
  drive_info <- shell('wmic logicaldisk get name', intern = TRUE)
  search_roots <- trimws(grep("^[A-Z]:", drive_info, value = TRUE))
  
  found_paths <- c()
  for (root_path in search_roots) {
    path_to_check <- file.path(root_path, folder_name)
    if (dir.exists(path_to_check)) {
      folders <- list.dirs(path_to_check, recursive = F, full.names = T)
      if(length(folders) == 0){
        found_paths <- c(found_paths, path_to_check)
      }else{
        found_paths <- c(found_paths, folders)  
      }
      
    }
  }
  return(found_paths)
}

## Three high Europe. Global no tiles. https://glad.umd.edu/users/Potapov/Europe_TCH/Tree_Height/. See _utilities_download_covariates.R code

tree_high_europe <- find_root_folder("tree_high_Europe")
full_path <- list.files(tree_high_europe, "*.tif$", full.names = T, recursive = F)

tree_high_europe <-
  tibble(path = full_path) |>
  tidyr::extract(
    col = path,
    into = c("var", "year"),
    regex = ".*/(.+?)_(\\d{4}).*\\.tif$",
    convert = TRUE,
    remove = FALSE
  ) |>
  dplyr::mutate(
    dataset = "tree_high_eu",
    UTM200 = NA,
    UTM10 = NA
  ) |>
  dplyr::select(dataset, var, path, year, UTM200, UTM10) |> 
  arrange(year)


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
  arrange(year)

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
  arrange(year)

master_covariates <- bind_rows(tree_high_europe, luke, latvus)
rm(tree_high_europe, luke, latvus)
#--------------------

# Prepare METSO and NON_METSO stands and routes polygons

coords <- st_read("data/coords_transformed.gpkg")
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

coords_utm_join <- coords |> 
  select(vakio) |> 
  mutate(set = "coords") |> 
  st_transform(st_crs(utm_200)) |> 
  st_join(y = utm_200) |> 
  rename(UTM200 = lehtitunnus) |> 
  st_join(y = utm_10) |> 
  rename(UTM10 = lehtitunnus) |> 
  select(vakio, set, UTM200, UTM10) |> 
  st_buffer(dist = 150)

# Add a unique ID to each polygon set to ensure results join correctly

metso_utm_join$poly_id <- 1:nrow(metso_utm_join)
coords_utm_join$poly_id <- 1:nrow(coords_utm_join)

rm(utm_10, utm_200)

#-----------------------

# Actually Get covariates information 

## Creates a direct-readable path for both zipped and unzipped rasters.
prepare_raster_paths <- function(df) {
  df |>
    mutate(
      path = gsub("\\\\", "/", path),
      proc_path = if_else(
        endsWith(path, ".zip"),
        paste0("/vsizip/", path, "/", sub("\\.zip$", "", basename(path))),
        path
      )
    ) 
}

## Gets the CRS from a raster file as a WKT string.

get_raster_crs <- function(raster_path) {
  tryCatch({
    r <- terra::rast(raster_path)
    return(terra::crs(r, proj = TRUE))
  }, error = function(e) {
    warning(paste("Could not read CRS from:", raster_path, "\nError:", e$message))
    return(NA_character_)
  })
}

master_covariates <- master_covariates |>
  prepare_raster_paths() |>
  group_by(dataset) |>
  mutate(raster_crs = get_raster_crs(first(proc_path))) |>
  ungroup() |>  
  filter(!is.na(raster_crs))

## Extracts raw information values from covarites for a polygon or set of polygons
## Iterates over a data.frame of covariates information

##raster_info_df
## dataset (name of the set, character)
## var (specific name of the variable, character)
## file (path to the raster layer, character)
## year (numeric)
## UTM200 (UTM 200 code if the raster is divided in tiles, character)
## UTM10 (UTM 10 code if the raster is divided in tiles, character)

## polygons_sf
## set (type of polygon, character)
## UTM200 (location in UTM 200 tiles of the polygon, character)
## UTM10 (location in UTM 10 tiles of the polygon, character)

## value (path of RDS containing raw information)

extract_raw_values <- function(polygons_sf, raster_info_df, folder = "results") {
  
  message(paste("Fetching from", raster_info_df$dataset[1], "variable", unique(raster_info_df$var), " with polygon set", unique(polygons_sf$set)))
  
  target_crs <- raster_info_df$raster_crs[1]
  
  polygons_reprojected <- st_transform(polygons_sf, crs = target_crs, quiet = TRUE)
  
  folder_name <- file.path(folder, "intermediate_by_tile", paste0(raster_info_df$dataset[1], "-", unique(polygons_sf$set)))
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  
  created_files <- c() 
  
  for (i in seq_len(nrow(raster_info_df))) {
    
    raster_info <- raster_info_df[i, ]
    
    if (is.na(raster_info$UTM200) && is.na(raster_info$UTM10)) {
      polygons_to_process <- polygons_reprojected
      nm <- paste0(raster_info$var, "_", raster_info$year, "_", i, "_tile.rds")
    } else if (!is.na(raster_info$UTM200)) {
      polygons_to_process <- polygons_reprojected |> filter(UTM200 == raster_info$UTM200)
      nm <- paste0(raster_info$var, "_", raster_info$year, "_", raster_info$UTM200, "_tile.rds")
    } else if (!is.na(raster_info$UTM10)) {
      polygons_to_process <- polygons_reprojected |> filter(UTM10 == raster_info$UTM10)
      nm <- paste0(raster_info$var, "_", raster_info$UTM10, "_tile.rds")
    } else {
      next # Skip to the next iteration
    }
    
    if (nrow(polygons_to_process) == 0) {
      next # Skip to the next iteration
    }
    
    suppressWarnings({
      raw_data_for_tile <- exact_extract(
        x = rast(raster_info$proc_path),
        y = polygons_to_process,
        fun = NULL,
        include_cols = "poly_id"
      )
    })
    
    raw_data_for_tile <- bind_rows(raw_data_for_tile) |>
      mutate(variable = raster_info$var,
             dataset = raster_info$dataset,
             year = raster_info$year)
    
    tile_filename <- file.path(folder_name, nm)
    
    saveRDS(raw_data_for_tile, file = tile_filename)
    
    created_files <- c(created_files, tile_filename)
    
    rm(raw_data_for_tile); gc()
  }
  
  return(created_files)
}

list_of_groups <- master_covariates |> 
  group_by(dataset, var) |> 
  group_split()

coords_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(coords_utm_join, .x))
coords_raw_data <- purrr::map_dfr(compact(coords_raw_data), readRDS)

metso_raw_data <- purrr::map(list_of_groups, ~extract_raw_values(metso_utm_join, .x))
metso_raw_data <- purrr::map_dfr(compact(coords_raw_data), readRDS)