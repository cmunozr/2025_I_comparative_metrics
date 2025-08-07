library(dplyr)
library(sf)
library(terra)

# shape de corte
path_shape <- "data/finland_square.gpkg"

mask_ <- path_shape %>% 
  sf::st_read() %>%  
  sf::st_transform(st_crs(4326)) %>% 
  terra::vect()

dirs <- "D:/tree_high_Europe/"
new_dir <- paste0(dirname(dirs), basename(dirs), "_fi")
dir.create(new_dir, showWarnings = F, recursive = T)

for(i in 1:length(dirs)){
  #i <- 1
  message(dirs[i])
  listi <- list.files(dirs, recursive = F, full.names = T, pattern = "*.adf$|*.tif$")
  
  for(a in 1:length(listi)){
    #a <- 1
    message(listi[a])
    ras.a <- terra::rast(listi[a])
    print(listi[a])
    mask.a <- ras.a %>% terra::crop(mask_)
    terra::writeRaster(
      x = mask.a, filename = paste0(new_dir[i], "/", basename(listi[a])),
      overwrite = T, NAflag = -9999, gdal = c("COMPRESS=LZW", "TFW=NO")
    )
  }
  
  rm(mask.a)
  gc()
  # dir_raster <- rasterOptions()[["tmpdir"]]
  # file.remove(list.files(dir_raster, full.names = T, recursive = F))
}

