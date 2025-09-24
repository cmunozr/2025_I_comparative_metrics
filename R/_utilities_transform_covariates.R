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

#-----------------

prepare_covariates_data_toplot <- function(covar_data) {
  
  # Add category labels for plotting
  plot_data <- covar_data |>
    mutate(category = case_when(
        polygon_source == "coords" ~ "route_buffer",
        metso == 1                   ~ "metso",
        metso == 0                   ~ "no_metso"
        )
      ) |> 
    filter(!is.na(category)) |>
    mutate(category = factor(category, levels = c("route_buffer", "metso", "no_metso")))
  
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
  
  output_filename <- file.path(folder, paste0("plot_", covariate_name, ".jpeg"))
  ggsave(output_filename, p, width = 12, height = 8, dpi = 300)
  message(paste("Saved plot to:", output_filename))
  
}

#-----------------

generate_binary_plot <- function(covariate_name, data, folder) {
  message(paste("Generating BAR CHART for binary covariate:", covariate_name))
  
  plot_data <- data |>
    filter(variable == covariate_name) |>
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
  
  output_filename <- file.path(folder, paste0("plot_", covariate_name, ".jpeg"))
  ggsave(output_filename, p, width = 12, height = 8, dpi = 300)
  message(paste("Saved plot to:", output_filename))
}

#-----------------
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

#' Find a data folder by searching across multiple hard disks/mount points.
#' This function works on Windows, Linux, and macOS.
#'
#' @param folder_name The name of the directory to search for (e.g., "tree_high_Europe").
#' @return The full path to the folder if found, otherwise it throws an error.

find_root_folder <- function(folder_name, opt = output_base_path )  {
  
  os <- Sys.info()['sysname']
  
  found_path <- NULL
  
  if (os == "Windows") {
    drive_info <- shell('wmic logicaldisk get name', intern = TRUE)
    search_roots <- trimws(grep("^[A-Z]:", drive_info, value = TRUE))
    
    for (root_path in search_roots) {
      path_to_check <- file.path(root_path, folder_name)
      if (dir.exists(path_to_check)) {
        found_path <- path_to_check
        break
      }
    }
    
  } else if (os %in% c("Linux", "Darwin")) {

    search_locations <- c(opt, "/mnt", "/media", "/srv", "~")
    
    potential_parent_dirs <- unlist(lapply(search_locations, function(loc) {
      suppressWarnings(list.dirs(loc, recursive = FALSE, full.names = TRUE))
    }))
    
    all_search_paths <- c(potential_parent_dirs, search_locations)
    
    for (path in all_search_paths) {
      full_path_to_check <- file.path(path, folder_name)
      if (dir.exists(full_path_to_check)) {
        found_path <- full_path_to_check
        break # Exit the loop once the folder is found
      }
    }
    
  } else {
    stop("Unsupported operating system: ", os)
  }
  
  if (!is.null(found_path)) {
    message(paste("Success! Found data at:", found_path))
    return(found_path)
  } else {
    stop(paste("Could not find the folder '", folder_name, "'. Please check the path.", sep=""))
  }
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


#-------------

#' Extract Raster Values for Polygons by Tile
#'
#' @description
#' Extracts pixel values from a set of raster files for a given set of polygons.
#' The process is optimized for memory efficiency by iterating through raster tiles
#' and saving intermediate results to disk.
#'
#' @details
#' This function matches polygons to their corresponding raster tiles (based on UTM
#' grid codes) before extraction. It includes a "dry run" mode via the `only_paths`
#' parameter, which allows you to generate the list of expected output file paths
#' without performing the computationally expensive extraction. This is useful for
#' planning and debugging workflows.
#'
#' @param polygons_sf An `sf` object containing the polygons for which to extract
#'   data. Must include columns for `year` and UTM tile identifiers (`UTM200`, `UTM10`).
#' @param raster_info_df A data frame with metadata for each raster file to be
#'   processed. It must include columns like `dataset`, `var`, `year`, `proc_path`
#'   (the file path to the raster), `raster_crs`, and UTM tile identifiers.
#' @param root_output_dir A string specifying the root directory where the output files
#'   (in an "intermediate_by_tile" subdirectory) will be saved.
#' @param id_column A string with the name of the column in `polygons_sf` that
#'   contains the unique polygon identifier to be included in the final output.
#' @param only_paths A logical value. If `TRUE`, the function will perform a "dry run"
#'   and return the expected output file paths without processing any data. If `FALSE`
#'   (the default), it will perform the full extraction and save the files.
#'
#' @return A character vector of the file paths that were created (if `only_paths = FALSE`)
#'   or that would be created (if `only_paths = TRUE`).
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming 'master_covariates' and 'coords_utm_join' are pre-loaded
#'   list_of_groups <- master_covariates |>
#'     dplyr::group_by(dataset, var, year) |>
#'     dplyr::group_split()
#'
#'   # Perform a dry run to get file paths without processing
#'   expected_files <- extract_raw_values(
#'     polygons_sf = coords_utm_join,
#'     raster_info_df = list_of_groups[[1]],
#'     only_paths = TRUE
#'   )
#'
#'   # Perform the full extraction
#'   created_files <- extract_raw_values(
#'     polygons_sf = coords_utm_join,
#'     raster_info_df = list_of_groups[[1]],
#'     only_paths = FALSE
#'   )
#' }

extract_raw_values <- function(polygons_sf, 
                               raster_info_df, 
                               root_output_dir, 
                               id_column, 
                               only_paths) {
  
  message(paste("Fetching from", raster_info_df$dataset[1], "variable", raster_info_df$var[1], "year",
                raster_info_df$year[1], "with polygon set", polygons_sf$set[1]))
  
  target_crs <- raster_info_df$raster_crs[1]
  
  if(polygons_sf$set[1] == "metso"){
    polygons_sf <- polygons_sf |> 
      mutate(year = raster_info_df$year[1])
    progress_bar <- T
  }else{
    progress_bar <- F
  }
  
  folder_name <- file.path(root_output_dir, "intermediate_by_tile", paste0(raster_info_df$dataset[1], "-", polygons_sf$set[1]))
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  
  if(!only_paths){
    polygons_reprojected <- polygons_sf |>
      filter(year == raster_info_df$year[1]) |>
      st_transform(crs = target_crs) 
  }else{
    polygons_reprojected <- polygons_sf |> 
      st_drop_geometry() |> 
      filter(year == raster_info_df$year[1]) 
  }
  
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
    
    tile_filename <- file.path(folder_name, nm)
    
    if(!only_paths){
      
      suppressWarnings({
        raw_data_for_tile <- exact_extract(
          x = rast(raster_info$proc_path),
          y = polygons_for_this_tile,
          fun = NULL,
          include_cols = "poly_id",
          progress = progress_bar
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
      
      saveRDS(raw_data_for_tile, file = tile_filename)
      rm(raw_data_for_tile);gc
      
    }

    created_files <- c(created_files, tile_filename)
    
    rm(polygons_for_this_tile); gc()
  }
  
  return(created_files)
}

#---------------------

# Check and split if the utm code actually has multiple utm 200 zones, example: KL2 
split_utm_code <- function(code) {
  if (str_detect(code, "[A-Z]{2,}")) {
    letters_part <- str_extract(code, "[A-Za-z]+")
    number_part <- str_extract(code, "\\d+")
    individual_letters <- str_split_1(letters_part, "")
    return(paste0(individual_letters, number_part))
  } else {
    return(code)
  }
}

#----------------
#' Aggregate and Save Covariate Data
#'
#' This function reads raw covariate data from a list of file paths, processes it,
#' standardizes variable names using a dictionary, and saves the final data frame
#' as an RDS file.
#'
#' @param raw_data_list A list where each element is a character vector of
#'   paths to RDS files for a specific covariate group.
#' @param polygon_type A character string, either "coords" or "metso",
#'   indicating the type of polygon data being processed. This is also used for
#'   naming the output folder and files.
#' @param dict_covar A data frame used as a dictionary to standardize
#'   covariate variable names. It should contain columns like 'file_name' and
#'   'var_shultz'.
#' @param output_root_dir A character string specifying the root directory
#'   where the aggregated files will be saved.
#' @param path_pattern A character string with the pattern in the file path to
#'   be replaced. Default is "results".
#' @param path_replacement A character string for the replacement in the file
#'   path. Default is "D:".
#'
#' @return The function is used for its side effect of writing files to disk and
#'   will invisibly return a character vector of the paths to the saved files.
#'
aggregate_covariates <- function(raw_data_list,
                                 polygon_type,
                                 dict_covar,
                                 output_root_dir,
                                 df) {
  
  folder_name <- file.path(output_root_dir, "intermediate_aggregated", polygon_type)
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  
  id_cols <- c("standid", "sampleUnit")
  saved_files <- c()
  
  for (i in 1:length(raw_data_list)) {
    message(paste("Processing", polygon_type, "group", i, "of", length(raw_data_list)))
    
    file_paths <- raw_data_list[[i]]
    
    dfa <- purrr::map_dfr(
      .x = file_paths,
      .f = ~ purrr::compact(readRDS(.x))
    )
    
    if (nrow(dfa) == 0) {
      message(paste("  - Skipping group", i, "as it contains no data."))
      next
    }
    
    dfa <- dfa |>
      # Remove rows with NA values or a specific GIS NA value (32766).
      dplyr::filter(!is.na(value), value != 32766) |>
      dplyr::rename_with(~"polygon_id", .cols = dplyr::any_of(id_cols))
    
    dfa_var <- dfa$variable[1]
    
    if (!is.na(dfa_var)) {
      if (dfa$dataset[1] == "luke") {
        index <- which(paste0(dfa_var, "_") == dict_covar$file_name)
        dfa$variable <- dict_covar$var_shultz[index]
      } else {
        index <- which(paste0(dfa_var) == dict_covar$file_name)
        dfa$variable <- dict_covar$var_shultz[index]
      }
    } else {
      dfa$variable <- "tree_removal_latest_date"
    }
    
    dfa$polygon_source <- polygon_type
    
    join_column_name <- if (polygon_type == "metso") "standid" else "sampleUnit"
    
    dfa <- semi_join(dfa, df, by = c("polygon_id" = join_column_name)) |>
      dplyr::mutate(polygon_id = as.factor(polygon_id))
    
    output_file_name <- paste0(
      dfa$variable[1], "_",
      dfa$year[1], "_",
      polygon_type,
      ".rds"
    )
    
    output_file_path <- file.path(folder_name, output_file_name)
    
    saveRDS(object = dfa, file = output_file_path)
    saved_files <- c(saved_files, output_file_path)
    
    rm(dfa)
    gc()
  }
  
  return(invisible(saved_files))
}

#----------------------
# Weighted standard deviation using 'coverage_fraction' as weight
# This function now expects 'x' (value) and 'w' (coverage_fraction)
weighted_sd <- function(x, w, na.rm = TRUE) {
  if (na.rm) {
    complete_cases <- !is.na(x) & !is.na(w)
    x <- x[complete_cases]
    w <- w[complete_cases]
  }
  if (length(x) == 0) return(NA) # Return NA if no valid data
  w_mean <- stats::weighted.mean(x, w, na.rm = FALSE)
  sqrt(sum(w * (x - w_mean)^2) / sum(w))
}

#----------------------
# Ratio calculation
# This function expects 'cov_frac' and a 'val' column (which is the 'value' column for these vars)
# The original logic returns the sum of coverage fractions where value is 1
ratio <- function(cov_frac, val, na.rm = TRUE) {
  if (na.rm) {
    complete_cases <- !is.na(cov_frac) & !is.na(val)
    cov_frac <- cov_frac[complete_cases]
    val <- val[complete_cases]
  }
  if (length(cov_frac) == 0) return(NA)
  sum(cov_frac[val == 1])
}

#------------------------
# Weighted mean using 'coverage_fraction' as weight
weighted_mean <- function(x, w, na.rm = TRUE) {
  if (na.rm) {
    complete_cases <- !is.na(x) & !is.na(w)
    x <- x[complete_cases]
    w <- w[complete_cases]
  }
  if (length(x) == 0) return(NA)
  stats::weighted.mean(x, w, na.rm = FALSE)
}

#----------------------
# Maximum value
max_ <- function(x, na.rm = TRUE) {
  if (length(x[!is.na(x)]) == 0) return(NA)
  max(x, na.rm = na.rm)
}

#-----------------------
# Binary indicator
bin <- function(x, aux, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) == 0) return(NA)
  # Check if ANY value meets the condition
  ifelse(any(x >= aux), 1, 0)
}

#----------------------
#' Determine the Old-Growth Age Threshold Based on Biotope Class
#'
#' This function extracts the primary biotope class from a polygon identifier
#' and returns the corresponding age threshold for old-growth forest classification.
#' The age thresholds are hard-coded based on established values for Finnish forest
#' types (Forsius et al., 2021; Shultz et al., 2025).
#'
#' @param polygon_id A character string representing the polygon's unique ID,
#'   where the last character is the biotope class code (e.g., "43_2009_K").
#'
#' @return A numeric value representing the age threshold in years.
#'
#' @details
#' The function uses a `case_when` statement to apply specific rules:
#'
#' - **K - KUUSIMETSÄ (Spruce Forest):** This category represents forests where
#'   spruce is the dominant tree species. The model applies a **100-year** age
#'   threshold, a direct classification that correctly identifies mature spruce
#'   stands based on established ecological standards for this forest type.
#'
#' - **M - MÄNTYMETSÄ (Pine Forest):** This code is for forests dominated by pine.
#'   The model applies a **120-year** age threshold. This classification ensures
#'   that only pine forests that have reached a significant age, characteristic
#'   of old-growth conditions for this species, are identified.
#'
#' - **S - SEKAMETSÄ (Mixed Forest):** This category is for mixed forests with a
#'   combination of coniferous and deciduous trees. To prevent overestimation in
#'   these ecologically complex habitats, a conservative approach is taken by applying
#'   the highest age threshold (**120 years**, equivalent to pine). This ensures
#'   only exceptionally old and structurally developed mixed stands are flagged.
#'
#' - **L & J - LEHTIMETSÄ & JALOPUUMETSÄ (Deciduous & Noble Broadleaf Forest):**
#'   These codes represent forests dominated by broadleaf species. They are grouped
#'   and assigned a **60-year** age threshold, correctly treating deciduous
#'   habitats as a distinct category with unique aging dynamics.
#'
#' - **P - PENSASTO (Shrubland):** This category includes areas dominated by shrubs
#'   and young trees. While a woody habitat, it generally lacks mature trees.
#'   Assigning the high **120-year "pine" threshold** is a conservative method
#'   that correctly results in a non-old-growth classification in almost all cases.
#'
#' - **R - RÄME, NEVA, LETTO (Pine Mire, Bog, Rich Fen):** This category represents
#'   peatland habitats often dominated by slow-growing, stunted pines. Applying
#'   the **120-year "pine" threshold** is an ecologically sound choice, accounting
#'   for the fact that trees in these mire systems must reach an advanced age to be
#'   considered part of an old-growth ecosystem.
#'
#' - **Other Biotopes:** Any other biotope code (e.g., `H` for clear-cut, `A` for
#'   agricultural land) is assigned an effectively impossible threshold of `9999`
#'   years, ensuring they are correctly classified as non-forest for this analysis.
#'
get_age_threshold <- function(polygon_id) {
  biotope_class <- stringr::str_sub(polygon_id, -1)
  case_when(
    biotope_class %in% c("M", "S", "P", "R") ~ 120,
    biotope_class == "K"                      ~ 100,
    biotope_class %in% c("L", "J")           ~ 60,
    TRUE                                      ~ 9999
  )
}

#----------------------
# Returns 1 if ANY pixel in the polygon meets its specific age threshold.
old_grow <- function(stand_age, polygon_id) {
  age_threshold <- get_age_threshold(polygon_id[1]) # All IDs in a group are the same
  is_old <- stand_age >= age_threshold
  return(as.integer(any(is_old, na.rm = TRUE)))
}

#----------------------
# Returns the proportion (0.0 to 1.0) of pixels that meet their age threshold.
old_grow_ratio <- function(stand_age, polygon_id) {
  age_threshold <- get_age_threshold(polygon_id[1])
  is_old <- stand_age >= age_threshold
  return(mean(is_old, na.rm = TRUE))
}

#----------------------
# Returns the sum of the ages of ONLY the pixels that meet their threshold.
old_grow_agesum <- function(stand_age, polygon_id) {
  age_threshold <- get_age_threshold(polygon_id[1])
  # Subset the ages that meet the condition before summing
  ages_to_sum <- stand_age[stand_age >= age_threshold]
  return(sum(ages_to_sum, na.rm = TRUE))
}