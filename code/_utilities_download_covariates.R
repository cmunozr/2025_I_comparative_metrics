library(rvest)
library(dplyr)
library(stringr)
library(xml2)
library(fs)
library(terra)
library(tools)
options(timeout = 1000)

dict_covar <- read.csv(file.path("data", "covariates", "dictionary_covariates.csv"), sep = ";")

years <- c(2009, 2011, 2013, 2015, 2017, 2019, 2021)

#-----------------
# Download TREE HEIGHT
# https://www.sciencedirect.com/science/article/pii/S0034425723003486?via%3Dihub

target_disk <- "D:"

base_url <- "https://glad.umd.edu/users/Potapov/Europe_TCH/Tree_Height/"

webpage <- read_html(base_url)

links <- webpage |> 
  html_nodes("a") |>
  html_attr("href")

tif_links <- links |> 
  as_tibble() |> 
  filter(str_detect(value, "\\.tif$")) |> 
  mutate(full_url = paste0(base_url, value)) |> # Construct full URL
  pull(full_url)

download_dir <- paste0(target_disk, "\tree_high_Europe")

dir.create(download_dir, showWarnings = F)

paths_in <- character()

for (url in tif_links) {
  # Extract the filename from the URL
  filename <- basename(url)
  # Define the full path for saving the file
  dest_file <- file.path(download_dir, filename)
  
  # Define the full path for saving the file
  download_path <- file.path(download_dir, filename)
  paths_in[url] <- download_path
  tryCatch({
    download.file(url, destfile = dest_file, mode = "wb") # 'wb' for binary files
  }, error = function(e) {
    message(paste0("Error downloading ", filename, ": ", e$message))
  })
}

#-----------------
# Latvus model is quite big, it is needed a way to identify space disk in my storage units

get_all_disk_space_windows <- function() {
  command <- "wmic"
  args <- c("logicaldisk", "get", "Caption,Freespace,Size")
  result <- system2(command, args, stdout = TRUE, stderr = TRUE)
  
  # Filter out empty lines, headers, and any "No Instance(s)" messages
  result <- result[nchar(result) > 0]
  result <- result[!grepl("Caption", result) & !grepl("---", result) & !grepl("No Instance", result)]
  
  parsed_lines <- lapply(result, function(line) {
    parts <- str_split(line, "\\s\\s+")[[1]]
    parts <- parts[parts != ""]
    return(parts)
  })
  
  df_output <- do.call(rbind, parsed_lines) |> as.data.frame(stringsAsFactors = FALSE)
  colnames(df_output) <- c("Caption", "FreeSpace_bytes", "Size_bytes")
  
  df_output$FreeSpace_bytes <- as.numeric(df_output$FreeSpace_bytes)
  df_output$Size_bytes <- as.numeric(df_output$Size_bytes)
  
  # Filter out any rows that might have non-numeric data or a carriage return after conversion
  df_output <- df_output |>
    filter(!is.na(FreeSpace_bytes) & !is.na(Size_bytes)) |>
    filter(!if_any(everything(), ~ str_detect(., "\r")))
  
  return(df_output)
}

#---------------------
# Download from metsakeskus: LATVUS model

preferred_disk_order <- c("D:", "E:")

base_url <- "https://avoin.metsakeskus.fi/aineistot/Latvusmalli/Karttalehti/"

webpage <- read_html(base_url)

links <- webpage |> 
  html_nodes("a") |>
  html_attr("href")

links <- links[str_which(links, pattern = "\\d")]

folder_links <- links |> 
  as_tibble() |> 
  mutate(full_url = paste0(base_url, value)) |> # Construct full URL
  pull(full_url)


for(i in 1:length(folder_links)){
  #i <- 1
  link <- folder_links[i]
    
  subwebpage <- read_html(link)
  
  # first, lets to check if  we have enough storage space to download
  
  ## How many bites do we need?
  pre_content <- subwebpage |>
    html_node("pre") |> 
    html_text() 
  
  lines <- str_split(pre_content, pattern = " ")[[1]]
  
  index_size <- lines |> 
    str_which("CHM") -1
  
  total_size <- lines[index_size] |> 
    as.numeric() |> 
    sum()
  total_size <- total_size + (10000*1024)
  
  #--
  
  all_disk_info <- get_all_disk_space_windows()
  
  # Prioritize by preferred_disk_order, then by most free space
  candidate_disks <- all_disk_info |>
    mutate(
      preferred_order = match(Caption, preferred_disk_order), # Assign order based on preference
      preferred_order = ifelse(is.na(preferred_order), nrow(all_disk_info) + 1, preferred_order) # Put non-preferred at the end
    ) |>
    arrange(preferred_order, desc(FreeSpace_bytes)) # Sort by preference, then by most free space
  
  for (k in 1:nrow(candidate_disks)) {
    disk_candidate <- candidate_disks$Caption[k]
    free_space_candidate <- as.numeric(candidate_disks$FreeSpace_bytes[k])
    
    if (free_space_candidate > total_size) {
      chosen_target_disk <- disk_candidate
      break # Found a suitable disk, exit the inner loop
    }
  }
  
  # If no suitable disk was found after checking all
  if (is.null(chosen_target_disk)) {
    message(paste0("No disk found with enough free space for ", total_size_needed_with_buffer, " bytes for folder ", basename(link), ". Skipping download.\n"))
    stop() # Skip to the next folder in folder_links
  } else {
    message(paste0("Enough space found on ", chosen_target_disk, ". Proceeding with download."))
  }

  yr_i <- paste(str_extract_all(links[i], "\\d")[[1]], collapse = "") |> 
    as.numeric()
  
  if(sum(yr_i %in% years) != 1){
    message(paste0("Year ", yr_i, " is not in the list of target years ", paste(years, collapse = " "), ". Skipping download.\n"))
    next
  }
    
  base_url_i <- paste0(base_url, yr_i, "/")
    
  links_subwebpage <- subwebpage |> 
    html_nodes("a") |>
    html_attr("href")
    
  tif_links_subwebpage <- links_subwebpage |> 
    as_tibble() |> 
    filter(str_detect(value, "\\.tif$")) |> 
    mutate(full_url = paste0(base_url_i, value)) |> # Construct full URL
    pull(full_url)
    
    
  download_dir <- file.path(chosen_target_disk, "latvus", yr_i)
    
  dir.create(download_dir, showWarnings = F, recursive = T)

  for (a in 1:length(tif_links_subwebpage)) {
    #a <- 1
    url <- tif_links_subwebpage[a]
    # Extract the filename from the URL
    filename <- basename(url)
    # Define the full path for saving the file
    dest_file <- file.path(download_dir, filename)
    
    if (file.exists(dest_file)) {
      message(paste0("File already exists: ", filename, ". Skipping download."))
      next # Skip download if file exists
    }
    
    tryCatch({
      download.file(url, destfile = dest_file, mode = "wb") # 'wb' for binary files
    }, error = function(e) {
      message(paste0("Error downloading ", filename, ": ", e$message))
    })
  }
  message(paste0("\n"))
}

#---------------------
# download from LUKE: 
# http://www.nic.funet.fi/index/geodata/luke/vmi/
# Average stand diameter, Average stand length, Canopy cover broadleaves
# Canopy cover whole stand, Stand age, Stand basal area, Total wood volume
# Spruce volume, pine volume, birch volume, broadleaves volume

dict_luke <- dict_covar |> 
  filter(set == "luke", to_download == 1, downloaded == 0)

base_url <- "http://www.nic.funet.fi/index/geodata/luke/vmi/"

prefered_disk <- "D:"

for(i in 1:nrow(dict_luke)){
  # i <- 1
  var_name <- dict_luke[i, "file_name"]
  
  url_years <- paste0(base_url, years)
  
  for(a in 1:length(url_years)){
    # a <- 2
    webpage <- read_html(url_years[a])
    
    links <- webpage |> 
      html_nodes("a") |>
      html_attr("href") |> 
      str_subset(".tif")
    
    var_name_ <- str_replace(string = var_name, pattern = "(.*)_([^_]*)$", replacement = "\\1\\2")
    
    link_detected <- links[str_detect(links, pattern = paste0("^", var_name_, "(\\.|_)"))]
    
    link_to_download <- paste0(url_years[a], "/", link_detected)
    
    download_dir <- file.path(prefered_disk, "luke_", years[a])
    
    dir.create(download_dir, showWarnings = F, recursive = T)
    
    # Define the full path for saving the file
    dest_file <- file.path(download_dir, paste0(var_name_, ".tif"))
    
    tryCatch({
      download.file(link_to_download, destfile = dest_file, mode = "wb") # 'wb' for binary files
    }, error = function(e) {
      message(paste0("Error downloading ", filename, ": ", e$message))
    })
    
  }
  
}

#---------------------
# download from Maanmittauslaitos 
# https://www.nic.funet.fi/index/geodata/mml/dem10m/2019/
# Elevation model 10m

# Download TIF files recursively from a web directory
# Scans a starting URL, navigates through all subdirectories,
# and downloads any file ending with ".tif".
# url The URL of the directory to scan.
# output_dir The local folder where files will be saved.

download_tifs_recursive <- function(url, output_dir) {
  
  cat("Scanning directory:", url, "\n")
  
  tryCatch({
    webpage <- read_html(url)
    
    links <- webpage |>
      html_nodes("a") |>
      html_attr("href")
    
    for (link in links) {
      full_url <- paste0(url, link)
      
      if (str_ends(link, "/") && !str_starts(link, "/")) {
        # If the link is a directory (ends with "/"), call this function again
        # This is a "recursive" step!
        download_tifs_recursive(full_url, output_dir)
        
      } else if (str_ends(link, "\\.tif$")) {

        file_name <- basename(full_url)
        destination_path <- file.path(output_dir, file_name)
        
        # Check if the file already exists to avoid re-downloading
        if (!file.exists(destination_path)) {
          cat("-> Downloading:", file_name, "\n")
          download.file(full_url, destfile = destination_path, mode = "wb", quiet = TRUE)
        } else {
          cat("-> Skipping:", file_name, "(already exists)\n")
        }
      }
    }
    
  }, error = function(e) {
    # If a URL is not accessible, print an error message and continue
    cat("Error accessing", url, ":", e$message, "\n")
  })
}

start_url <- "https://www.nic.funet.fi/index/geodata/mml/dem10m/2019/"

output_folder <- "C:/DEM_Finland"
dir.create(output_folder)

download_tifs_recursive(url = start_url, output_dir = output_folder)

# making bigger tiles to speed process

dem_files <- list.files(output_folder, pattern = ".tif$", full.names = T, ) |> 
  tools::file_path_sans_ext() |> 
  basename() |> 
  str_extract_all("(\\D+\\d)|(\\d+)", simplify = T) |> 
  as.data.frame() |> 
  pull(V1) |> 
  data.frame("path" = dem_files, "utmzone" = _)

index <- unique(dem_files$utmzone)

for(i in 1:length(index)){
  tryCatch({
    message(paste("processing utm zone:", index[i]))
    dem_files |> 
      dplyr::filter(utmzone == index[i]) |> 
      pull(path) |> 
      sprc(lapply(terra::rast)) |> 
      merge(algo = 2) |> 
      terra::writeRaster(
        file.path(output_folder, paste0( index[i], ".tif")), 
        gdal = c("COMPRESS=DEFLATE"), 
        overwrite=TRUE
      )  
    }, error = function(e) {
      cat("Error accessing", index[i], ":", e$message, "\n")
    }
  )
}

file.remove(dem_files$path)


