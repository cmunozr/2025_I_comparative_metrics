library(sf)
library(dplyr)

paths_in <- list.files(path = "data/metso/raw2/", full.names = T, recursive = F, pattern = ".gpkg")
metso_trt_control <- read.csv(file.path("data", "metso", "treatment_control_standid.csv"))
epsg <-  "EPSG:3067"

gpkg <- read_sf(paths_in) |> 
  st_transform(crs = st_crs(epsg)) |> 
  st_zm() #2d only

metso_trt_control_stands <- gpkg[(gpkg$standid %in% metso_trt_control$standid), ] |> 
  left_join(metso_trt_control, by = "standid") |> 
  relocate(geom, .after = metso)

write_sf(metso_trt_control_stands, file.path("data", "metso", "treatment_control_stand_v3.gpkg"), delete_layer = T)




