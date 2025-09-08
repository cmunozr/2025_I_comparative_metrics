library(rvest)
library(dplyr)
library(stringr)
library(xml2)
library(fs)
options(timeout = 1000)

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

years <- c(2009, 2011, 2013, 2015, 2017, 2019, 2021)

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
# Spruce volume, pine volume, birch volume, 
# first run _utilities_download_shultz_covariates.py
# to activate environment in vscode: C:/.compMetrics/Scripts/Activate.ps1
# extract list of links from mail
# run next

luke_urls <- readClipboard() #comming from the mail
prefered_disk <- "E:"

for (i in 11:length(luke_urls)) {
  #i <- 10
  url_i <- luke_urls[i]
  
  webpage <- read_html(url_i)
  
  nm <- webpage |>
    html_node("table") |> 
    paste(concatenate = ",") |> 
    str_split(pattern = "<td>|</td>", simplify = T)
  nm <- (nm[2]) |> 
    str_replace_all(pattern = " |/|-", replacement = "_")

  yr_ <- str_extract_all(nm, pattern = "\\d")[[1]]
  
  if(length(yr_) == 5){
    yr_ <- yr_[-5]
  }
  
  yr_ <- yr_ |> 
    paste(collapse = "")
  
  links <- webpage |> 
    html_nodes("a") |>
    html_attr("href")
  
  download_dir <- file.path(prefered_disk, "luke", yr_)
  
  dir.create(download_dir, showWarnings = F, recursive = T)

  # Define the full path for saving the file
  dest_file <- file.path(download_dir, nm)
  
  dest_file_zip <- paste0(dest_file, "_batch1.zip")
  
  if (file.exists(dest_file_zip)) {
    dest_file_zip <- paste0(dest_file, "_batch2.zip")
  }
  
  tryCatch({
    download.file(links, destfile = dest_file_zip, mode = "wb") # 'wb' for binary files
  }, error = function(e) {
    message(paste0("Error downloading ", filename, ": ", e$message))
  })
}

#------

library(fs)

base_directory <- "E:/luke/"

# Get a list of all main folders (FolderX, e.g., Folder1, Folder2)
main_folders <- dir_ls(base_directory, type = "directory")

# Loop through each main folder
for (folder_path in main_folders) {
  #folder_path <- main_folders[1]
  folder_name <- path_file(folder_path)
  
  zip_files_in_folder <- dir_ls(folder_path, glob = "*.zip")
  
  for (zip_file_path in zip_files_in_folder) {
    # zip_file_path <- zip_files_in_folder[1]
    zip_file_name <- path_file(zip_file_path)
    
    # Create a temporary directory for extraction
    temp_extract_dir <- path(base_directory, paste0("temp_extract_", Sys.time() |> as.integer(), "_", sample(1:1000, 1)))
    dir_create(temp_extract_dir)
    
    # Attempt to unzip the current zip file
    unzip_success <- tryCatch({
      unzip(zip_file_path, exdir = temp_extract_dir, junkpaths = FALSE)
      TRUE
    }, error = function(e) {
      message(paste0("    ERROR unzipping '", zip_file_name, "': ", e$message, "\n"))
      FALSE
    })
    
    if (!unzip_success) {
      dir_delete(temp_extract_dir) # Clean up failed temp dir
      next
    }else{
      subfolder <- dir_ls(temp_extract_dir, type = "directory")
      
      
      # --- Process nested zip files ---
      nested_zip_files <- dir_ls(subfolder, glob = "*.zip")
      
      if (length(nested_zip_files) > 0) {
        for (nested_zip in nested_zip_files) {
          target_path <- path(folder_path, path_file(nested_zip))
          file_move(nested_zip, target_path)
        }
      }
      
      # --- Deleting ---
      
      file_delete(zip_file_path)
      dir_delete(temp_extract_dir)
    }
  }
}
