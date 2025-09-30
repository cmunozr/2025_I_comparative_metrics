library(sf)
library(dplyr)

paths_in <- list.files("data/metsakeskus/metsavarakuviot", full.names = T, recursive = F)
metso_trt_control <- read.csv("data/metso/treatment_control_standid.csv")

# Función para leer shapefiles desde archivos ZIP
leer_gpkg_desde_zip <- function(zip_file, epsg = "EPSG:3067", metso = metso_trt_control, tmp.fol = "data/metsakeskus/metsavarakuviot") {
  #zip_file <- paths_in[6]
  temp_dir <- tmp.fol # Gets a system-specific temporary directory
  gpkg_unzip_dir <- file.path(temp_dir, "unzipped_gpkg_data")
  
  unzip(zipfile = zip_file, exdir = gpkg_unzip_dir)
  gpkg_files <- list.files(gpkg_unzip_dir, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE)
  gpkg_dsn_unzipped <- gpkg_files[1]
  
  gpkg <- read_sf(gpkg_dsn_unzipped) |> 
    st_transform(crs = st_crs(epsg)) |> 
    st_zm() #2d only
  
  gpkg <- gpkg[(gpkg$standid %in% metso$standid), ] |> 
    left_join(metso, by = "standid") |> 
    relocate(metso, .before = geometry)
  
  unlink(gpkg_unzip_dir, recursive = TRUE)
  
  return(gpkg)
}

shps <- paths_in |> 
  lapply(X = ., function(X){
    tmp <- leer_gpkg_desde_zip(X)
  })

metso_trt_control_stands <- bind_rows(shps)

write_sf(metso_trt_control_stands, "data/metso/treatment_control_stand.gpkg", delete_layer = T)