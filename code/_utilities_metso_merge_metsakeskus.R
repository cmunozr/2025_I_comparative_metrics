library(sf)
library(dplyr)
library(stringr)

paths_in <- list.files(path = "data/metso/raw2/", full.names = T, recursive = F, pattern = ".gpkg")
metso_trt_control <- read.csv(file.path("data", "metso", "treatment_control_standid.csv")) |> 
  mutate(temp_standid = str_split(standid, pattern = "_", simplify = T)[,1])
epsg <-  "EPSG:3067"

gpkg <- read_sf(paths_in) |> 
  st_transform(crs = st_crs(epsg)) |> 
  st_zm() |> 
  mutate(standid = as.character(standid)) #2d only

metso_trt_control_stands <- gpkg[gpkg$standid %in% metso_trt_control$temp_standid, ] |> 
  left_join(metso_trt_control, join_by("standid" == "temp_standid") ) |> 
  relocate(geom, .after = metso) |>
  select(-standid) |>
  rename("standid" = "standid.y")

write_sf(metso_trt_control_stands, file.path("data", "metso", "treatment_control_stand_v4.gpkg"), delete_layer = T)




