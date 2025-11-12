library(sf)
library(ggplot2)
library(units)
library(stringr)
library(tidyr)
library(dplyr)
set.seed(11072024)

stands <- read_sf("data/metso/treatment_control_stand.gpkg")

distance_buffer_km <- 2

metso_stands <- stands |> 
  filter(metso == 1)

control_stands <- stands |> 
  filter(metso == 0) |> 
  sample_n(size = nrow(metso_stands), replace = F)

rm(stands); gc()

#-------- overlap 

coords <- read_sf("data/fbs/vakiolinja/Vakiolinjat_routes.geojson")

coords_buf <- coords |> 
  st_transform(crs = st_crs(metso_stands)) |>
  st_cast("POLYGON")|> 
  st_buffer(dist = distance_buffer_km*1000, endCapStyle = "FLAT") |> 
  mutate(vakio = str_split(name, pattern = ", ", simplify = T)[,1] |> 
           as.numeric())

over_coords_metso <- coords_buf[metso_stands, ] |> 
  sf::st_drop_geometry() |> 
  select(vakio) |> 
  mutate(stand = 1)

over_coords_control <- coords_buf[control_stands, ]|> 
  sf::st_drop_geometry() |> 
  select(vakio) |> 
  mutate(stand = 0)

#------- extract vakio id overlapping with stands

coords_id <- st_drop_geometry(coords_buf[, "vakio"])

metso_coords_id <- over_coords_metso |>
  distinct(vakio) |>
  mutate(is_metso = 1)

control_routes_id <- over_coords_control |>
  distinct(vakio) |>
  mutate(is_control = 1)

coords_dist_buf <- coords_id |>
  left_join(metso_coords_id, by = "vakio") |>
  left_join(control_routes_id, by = "vakio") |>
  mutate(
    is_metso = ifelse(is.na(is_metso), 0, is_metso),
    is_control = ifelse(is.na(is_control), 0, is_control)
  )

write.csv(coords_dist_buf, file.path("data", "fbs", paste0("route_distance_",distance_buffer_km, "km_metso_control.csv")), row.names = F)

#------- creating spatial object for visualization of the distribution

bind_cols(coords, coords_dist_buf) |> 
  sf::write_sf(file.path("data", "fbs", paste0("route_distance_", distance_buffer_km,"km_metso_control.gpkg")), delete_layer = T)
