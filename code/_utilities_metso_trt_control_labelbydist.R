# Some control stands are far away from metso sites, take them as it could be a
# bias positively the change

library(sf)
library(ggplot2)
library(units)
library(stringr)
library(tidyr)
library(dplyr)
set.seed(11072024)

stands <- read_sf(file.path("data", "metso", "treatment_control_stand_v3.gpkg"))

distance_buffer_km <- 2

metso_stands <- stands |> 
  filter(metso == 1)

control_stands <- stands |> 
  filter(metso == 0)

#-------- buffer metso and overlap control 

metso_stands_buf <- metso_stands |> 
  st_buffer(dist = distance_buffer_km*1000, endCapStyle = "ROUND") |> 
  summarise()
write_sf(metso_stands_buf, file.path("data", "metso", "metso_stands_buf.gpkg"))

control_stands_filtered <- control_stands[metso_stands_buf, ] |> 
  pull(standid)

stands <- stands |> 
  mutate( in_2km_dist  = ifelse(standid %in% control_stands_filtered, yes = 1, no = 0))

#------- creating spatial object

stands |> 
  sf::write_sf(file.path("data", "metso", paste0("treatment_control_stand_V3.gpkg")), delete_layer = T)


