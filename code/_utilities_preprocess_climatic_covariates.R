library(terra)
library(purrr)
library(tidyr)
library(zip)
library(sf)
library(stringr)
library(dplyr)

os <- Sys.info()['sysname']

if (os == "Windows") {
  output_base_path <- "D:"
} else if (os %in% c("Linux", "Darwin")) {
  output_base_path <- "/media/pitm/My Passport/"
}  

source(file.path("code","_utilities_transform_covariates.R"))

utm_200LR <- st_read(file.path("data", "utm35_zones", "TM35_karttalehtijako.gpkg"), layer = "utm200LR") |> 
  summarise(geometry = st_union(geometry), .groups = "drop")
  
clim <- find_root_folder("nordic_climate")
zip_path <- list.files(clim, full.names = T)

files_in_zip <- zip_list(zip_path) |> 
  dplyr::select(filename) |> 
  mutate(
    str_split(tools::file_path_sans_ext(filename), pattern = "_", simplify = T) |> 
      as.data.frame() |> 
      select(V2, V6)
  ) |> 
  rename(
    var = V2,
    date = V6
  ) |> 
  tidyr::extract(
    col = date,
    into = c("year", "month_day"),
    regex = "^(\\d{4})(\\d{4})$", # divide a string in two and arrange in two different columns
    convert = TRUE,
    remove = FALSE
  ) |> 
  group_by(var, year)

files_in_zip |>
  group_by(var, year) |>
  dplyr::group_walk(.f = ~ process_climate(..1, ..2, out_dir = clim, zip_path = zip_path))

