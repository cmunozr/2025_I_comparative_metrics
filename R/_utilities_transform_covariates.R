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

# library(dplyr)
# library(sf)
# library(terra)
# 
# # shape de corte
# path_shape <- "data/finland_square.gpkg"
# 
# mask_ <- path_shape %>% 
#   sf::st_read() %>%  
#   sf::st_transform(st_crs(4326)) %>% 
#   terra::vect()
# 
# dirs <- "D:/tree_high_Europe/"
# new_dir <- paste0(dirname(dirs), basename(dirs), "_fi")
# dir.create(new_dir, showWarnings = F, recursive = T)
# 
# for(i in 1:length(dirs)){
#   #i <- 1
#   message(dirs[i])
#   listi <- list.files(dirs, recursive = F, full.names = T, pattern = "*.adf$|*.tif$")
#   
#   for(a in 1:length(listi)){
#     #a <- 1
#     message(listi[a])
#     ras.a <- terra::rast(listi[a])
#     print(listi[a])
#     mask.a <- ras.a %>% terra::crop(mask_)
#     terra::writeRaster(
#       x = mask.a, filename = paste0(new_dir[i], "/", basename(listi[a])),
#       overwrite = T, NAflag = -9999, gdal = c("COMPRESS=LZW", "TFW=NO")
#     )
#   }
#   
#   rm(mask.a)
#   gc()
# }

