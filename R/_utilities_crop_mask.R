library(dplyr)
library(sf)
library(terra)

# shape de corte
path_shape <- "data/nordic_countries.gpkg"


mask_ <- path_shape %>% 
  sf::st_read() %>%  
  sf::st_transform(st_crs(4326)) %>% 
  terra::vect()
template <- rast(resolution = c(0.008333333, 0.008333333), crs = crs(mask_ ), ext(mask_), vals = 1)

dirs <- list.dirs("data/covariates/", recursive = F, full.names = T)[-1]
#dirs_nm <- list.dirs("neotropico/ALPHA/", recursive = F, full.names = F)[-2]

for(i in 1:length(dirs)){
  print(dirs[i])
  #i <- 1
  listi <- list.files(dirs, recursive = F, full.names = T, pattern = "*.adf$|*.tif$")
  #listi <- dirs
  
  for(a in 1:length(listi)){
    #a <- 1
    ras.a <- terra::rast(listi[a])
    print(listi[a])
    mask.a <- ras.a %>% raster::crop(template) # %>% raster::mask(mask.t)
    #mask.a <- trunc(mask.a)
    terra::writeRaster(
      x = mask.a, filename = paste0(dirs[i], "/", names(ras.a), "_.tif"),
      #x = mask.a, filename = paste0(dirs[i], "/", dirs_nm[a], "_.tif"),
      overwrite = T, NAflag = -9999, gdal = c("COMPRESS=LZW", "TFW=NO")
    )
  }
  
  rm(mask.a)
  gc()
  # dir_raster <- rasterOptions()[["tmpdir"]]
  # file.remove(list.files(dir_raster, full.names = T, recursive = F))
}
