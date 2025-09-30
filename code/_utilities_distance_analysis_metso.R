library(sf)
library(dplyr)
library(ggplot2)
library(units)

metso_trt_control_stands <- read_sf("data/metso/treatment_control_stand.gpkg") |> 
  group_by(metso) |> 
  sample_n(size = 46078 , replace = F) |> 
  ungroup()

coords <- read_sf("data/fbs/vakiolinja/Vakiolinjat_routes.geojson") |> 
  st_transform(crs = "EPSG:2393")

buffers <- coords |> 
  st_buffer(dist = 150, endCapStyle = "FLAT") |> 
  st_transform(crs = "EPSG:3067")

write_sf(buffers, "data/metso/routes_buffer_line.gpkg", delete_layer = T)

pols <- coords |>
  st_cast("POLYGON")|> 
  st_buffer(dist = 150, endCapStyle = "FLAT")|> 
  st_transform(crs = "EPSG:3067")
write_sf(pols, "data/metso/routes_polBufferLine.gpkg", delete_layer = T)

pols1km <- coords |>
  st_cast("POLYGON")|> 
  st_buffer(dist = 1000, endCapStyle = "FLAT")|> 
  st_transform(crs = "EPSG:3067")
write_sf(pols1km, "data/metso/routes_polBuffer1kmLine.gpkg", delete_layer = T)

coords <- coords |> 
  st_transform(crs = "EPSG:3067")

a <- metso_trt_control_stands[coords,]
write_sf(a, "data/metso/a_lines_touch.gpkg")
table(a$metso)

b <- metso_trt_control_stands[buffers,]
write_sf(b, "data/metso/b_bufferLine_touch.gpkg", delete_layer = T)
table(b$metso)

c <- metso_trt_control_stands[pols,]
write_sf(c, "data/metso/c_polBufferLine.gpkg", delete_layer = T)
table(c$metso)

d <- metso_trt_control_stands[pols1km,]
write_sf(d, "data/metso/d_polBuffer1kmLine.gpkg", delete_layer = T)
table(d$metso)

#---------

xy <- metso_trt_control_stands |>
  st_centroid()

nd <- st_nearest_feature(xy, buffers)

x <- buffers[nd, ]

distances <- st_distance(xy, x, by_element = TRUE)

nd <- data.frame("standid" = metso_trt_control_stands$standid, 
                 "distance" = distances, 
                 "metso" = as.factor(metso_trt_control_stands$metso))

nd |> 
  ggplot() +
  geom_density(aes(x = distance, fill = metso), alpha = 0.4) +
  theme_minimal()
