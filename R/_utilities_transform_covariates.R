#' Crop multiple rasters based on a vector mask
#'
#' This function scans a source directory for raster files, crops each one using
#' a provided vector file (e.g., shapefile, geopackage), and saves the cropped
#' output to a new directory.
#'
#' @param source_dir A string with the path to the directory containing the original raster files.
#' @param mask_path A string with the path to the vector file (e.g., .gpkg, .shp) to be used as a cropping mask.
#' @param output_dir A string with the path to the directory where cropped rasters will be saved. The directory will be created if it doesn't exist.
#' @param file_pattern A regular expression string to identify the raster files.
#'   Defaults to finding .tif and .adf files.
#'
#' @return Returns a character vector of the full paths to the newly created cropped raster files.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Define paths
#'   my_rasters_folder <- "C:/data/raw_rasters"
#'   my_mask_file <- "C:/data/study_area.shp"
#'   my_output_folder <- "C:/data/cropped_rasters"
#'
#'   # Run the function
#'   created_files <- crop_rasters_by_mask(
#'     source_dir = my_rasters_folder,
#'     mask_path = my_mask_file,
#'     output_dir = my_output_folder
#'   )
#'
#'   # See the list of new files
#'   print(created_files)
#' }
crop_rasters_by_mask <- function(source_dir,
                                 mask_path,
                                 output_dir,
                                 file_pattern = "\\.adf$|\\.tif$") {
  
  if (!requireNamespace("terra", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    stop("This function requires the 'terra' and 'sf' packages. Please install them.")
  }
  
  if (!dir.exists(source_dir)) {
    stop("Source directory not found: ", source_dir)
  }
  if (!file.exists(mask_path)) {
    stop("Mask file not found: ", mask_path)
  }
  
  mask_vect <- sf::st_read(mask_path, quiet = TRUE) |>
    terra::vect()
  
  # Create the output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find all raster files in the source directory that match the pattern
  raster_files <- list.files(
    source_dir,
    recursive = FALSE,
    full.names = TRUE,
    pattern = file_pattern
  )
  
  if (length(raster_files) == 0) {
    warning("No raster files found in the source directory matching the pattern.")
    return(invisible(NULL))
  }

  # This vector will store the paths of the new files
  output_files <- c()
  
  for (current_file in raster_files) {
    message(paste("Processing:", basename(current_file)))
    
    # Use a tryCatch block to handle potential errors with a file
    # and allow the loop to continue with the next file.
    tryCatch({
      original_raster <- terra::rast(current_file)
      
      # Reproject the mask to match the raster's CRS
      mask_reprojected <- terra::project(mask_vect, terra::crs(original_raster))
      
      cropped_raster <- terra::crop(original_raster, mask_reprojected, mask = TRUE)
      

      output_filename <- file.path(output_dir, basename(current_file))
      
      terra::writeRaster(
        x = cropped_raster,
        filename = output_filename,
        overwrite = TRUE,
        NAflag = -9999,
        gdal = c("COMPRESS=LZW", "TFW=NO") # GDAL options for compression
      )
      
      output_files <- c(output_files, output_filename)
      
    }, error = function(e) {
      # If an error occurs, print it and move on
      warning(paste("Failed to process", basename(current_file), ":\n", e$message))
    })
    
  }
  
  message("Processing complete!")
  return(output_files)
}

#-------------------
filter_covariate_data <- function(covariate_name, files, sampled_ids, route_map) {
  # Filter file paths to get only those relevant to the current covariate
  relevant_files <- files[str_detect(basename(files), covariate_name)]
  
  if (length(relevant_files) == 0) {
    message(paste("No data files found for covariate:", covariate_name))
    return(NULL)
  }
  
  # Map over the relevant files, reading and processing each one
  covariate_data_list <- lapply(relevant_files, function(file) {
    
    data_chunk <- readRDS(file)
    
    processed_chunk <- data_chunk |>
      filter(!is.na(value), value != 32766) |> # 32766 is NA in some GIS systems 
      rename_with(~"polygon_id", .cols = any_of(c("standid", "vakio")))
    
    if (!"polygon_id" %in% names(processed_chunk)) return(NULL)
    
    source_type <- processed_chunk$polygon_source[1]
    if (!is.na(source_type) && source_type == "metso") {
      final_chunk <- processed_chunk %>% filter(polygon_id %in% sampled_ids)
    } else if (!is.na(source_type) && source_type == "coords") {
      final_chunk <- semi_join(processed_chunk, route_map, by = c("polygon_id" = "vakio", "year"))
    } else {
      final_chunk <- NULL
    }
    return(final_chunk)
  })
  
  bind_rows(covariate_data_list)
}

#-----------------

prepare_covariates_data_toplot <- function(filtered_data, metso_map) {
  # Check for valid input
  if (is.null(filtered_data) || nrow(filtered_data) == 0) {
    return(NULL)
  }
  
  # Add category labels for plotting
  plot_data <- filtered_data %>%
    left_join(metso_map, by = c("polygon_id" = "standid")) %>%
    mutate(category = case_when(
      polygon_source == "coords" ~ "route buffer",
      metso == 1                   ~ "metso",
      metso == 0                   ~ "no metso"
    )) %>%
    filter(!is.na(category)) %>%
    mutate(category = factor(category, levels = c("route buffer", "metso", "no metso")))
  
  return(plot_data)
}

#-----------------

generate_covariate_plot <- function(covariate_name, data, folder) {

  message(paste("Generating violin-plot for:", covariate_name))
  
  plot_data <- data |>
    filter(variable == covariate_name, !is.na(value))
  
  # Ensure there's data to plot
  if (nrow(plot_data) == 0) {
    message(paste("No data available for", covariate_name, "- skipping."))
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = category, y = value, fill = category)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.3, outlier.shape = NA) +
    facet_wrap(~year, scales = "free_y") +
    labs(
      title = paste("Distribution of", covariate_name, "by Year"),
      subtitle = "Comparison across Route Buffers, METSO, and Non-METSO Stands",
      x = "Category",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  output_filename <- file.path(folder, paste0("plot_", covariate_name, ".png"))
  ggsave(output_filename, p, width = 12, height = 8, dpi = 300)
  message(paste("Saved plot to:", output_filename))
  
}

#-----------------

generate_binary_plot <- function(covariate_name, data, folder) {
  message(paste("Generating BAR CHART for binary covariate:", covariate_name))
  
  plot_data <- data %>%
    filter(variable == covariate_name) %>%
    mutate(value = as.factor(value))
  
  if (nrow(plot_data) == 0) {
    message("No data to plot.")
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = category, fill = value)) +
    geom_bar(position = "fill") +
    facet_wrap(~year) +
    labs(
      title = paste("Proportional Distribution of", covariate_name, "by Year"),
      subtitle = "Comparison across Route Buffers, METSO, and Non-METSO Stands",
      x = "Category",
      y = "Proportion",
      fill = "Value" # This will be the legend title
    ) +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percentages
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  output_filename <- file.path(folder, paste0("plot_", covariate_name, ".png"))
  ggsave(output_filename, p, width = 12, height = 8, dpi = 300)
  message(paste("Saved plot to:", output_filename))
}

#---------------------

# To calculate weighted standard deviation
weighted_sd <- function(x, w, na.rm = TRUE) {
  if (na.rm) {
    complete_cases <- !is.na(x) & !is.na(w)
    x <- x[complete_cases]
    w <- w[complete_cases]
  }
  # Calculate the weighted mean
  w_mean <- weighted.mean(x, w)
  # Calculate the weighted standard deviation
  sqrt(sum(w * (x - w_mean)^2) / sum(w))
}

#' Extends an sf LINESTRING by a specified distance.
#'
#' @param line An sf object containing a single LINESTRING.
#' @param distance The distance to extend the line by (in the line's coordinate system units).
#' @param end Which end to extend: "end" (default) or "start".
#' @return A new sf object with the extended LINESTRING.

#---------------------

st_extend_line <- function(line, distance, end = "end") {
  
  require(sf)
  require(dplyr)
  
  
  # Get the coordinates of the line's vertices
  xy <- st_coordinates(line)[, 1:2]
  
  # Determine the segment to use for calculating the angle
  if (end == "end") {
    # Use the last two vertices for the final segment
    segment_start <- xy[nrow(xy) - 1, 1:2]
    segment_end <- xy[nrow(xy), 1:2]
  } else {
    # Use the first two vertices for the initial segment
    segment_start <- xy[2, 1:2]
    segment_end <- xy[1, 1:2] # Reversed to get the outward angle
  }
  
  # Calculate the angle of the segment
  # atan2(y, x) gives the angle in radians
  angle <- atan2(segment_end[2] - segment_start[2], segment_end[1] - segment_start[1])
  
  # Calculate the coordinates of the new point using trigonometry
  new_point_x <- segment_end[1] + distance * cos(angle)
  new_point_y <- segment_end[2] + distance * sin(angle)
  new_point <- c(new_point_x, new_point_y)
  
  # Combine the new point with the original coordinates
  if (end == "end") {
    new_xy <- rbind(xy, new_point)
  } else {
    new_xy <- rbind(new_point, xy)
  }
  
  # Create the new, extended linestring
  extended_line <- st_linestring(new_xy) |>
    st_sfc(crs = st_crs(line)) |>
    st_sf() |> 
    cbind(st_drop_geometry(line))
  
  return(extended_line)
}

#--------------------------

#### Functions used to organize and fetch data from the raw covariates

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

#--------------

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

#-------------

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

extract_raw_values <- function(polygons_sf, raster_info_df, folder = "D:", id_column = "vakio") {
  
  message(paste("Fetching from", raster_info_df$dataset[1], "variable", raster_info_df$var[1], "year",
                raster_info_df$year[1], "with polygon set", polygons_sf$set[1]))
  
  target_crs <- raster_info_df$raster_crs[1]
  
  folder_name <- file.path(folder, "intermediate_by_tile", paste0(raster_info_df$dataset[1], "-", polygons_sf$set[1]))
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  
  polygons_reprojected <- polygons_sf |>
    filter(year == raster_info_df$year[1]) |>
    st_transform(crs = target_crs) 
  
  created_files <- c()
  
  for (i in 1:nrow(raster_info_df)) {
    
    message(paste("tile", i, "of", nrow(raster_info_df)))
    
    raster_info <- raster_info_df[i, ]
    
    if (!is.na(raster_info$UTM200)) {
      polygons_for_this_tile <- polygons_reprojected |> 
        filter(UTM200 == raster_info$UTM200)
      nm <- paste0(raster_info$var, "_", raster_info$year, "_", raster_info$UTM200, "_tile.rds")
      
    } else if (!is.na(raster_info$UTM10)) {
      polygons_for_this_tile <- polygons_reprojected |> 
        filter(UTM10 == raster_info$UTM10)
      nm <- paste0(raster_info$var, "_", raster_info$year, "_", raster_info$UTM10, "_tile.rds")
      
    } else {
      # This handles cases where the raster is not tiled (covers the whole area).
      polygons_for_this_tile <- polygons_reprojected
      nm <- paste0(raster_info$var, "_", raster_info$year, "_tile.rds")
      
    }
    
    # If no polygons fall within this tile, skip to the next raster.
    if (nrow(polygons_for_this_tile) == 0) {
      next
    }
    
    suppressWarnings({
      raw_data_for_tile <- exact_extract(
        x = rast(raster_info$proc_path),
        y = polygons_for_this_tile,
        fun = NULL,
        include_cols = "poly_id",
        progress = F
      )
    })
    
    raw_data_for_tile <- bind_rows(raw_data_for_tile) |>
      mutate(variable = raster_info$var,
             dataset = raster_info$dataset,
             year = raster_info$year) |>
      left_join(
        st_drop_geometry(polygons_reprojected) |> select(all_of(id_column), poly_id),
        by = "poly_id"
      )
    
    
    tile_filename <- file.path(folder_name, nm)
    
    saveRDS(raw_data_for_tile, file = tile_filename)
    
    created_files <- c(created_files, tile_filename)
    
    rm(raw_data_for_tile, polygons_for_this_tile); gc()
  }
  
  return(created_files)
}
