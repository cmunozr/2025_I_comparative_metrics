library(rvest)
library(dplyr)
library(stringr)
library(xml2)
options(timeout = 600)

base_url <- "https://avoin.metsakeskus.fi/aineistot/Metsavarakuviot/Maakunta/"

webpage <- read_html(base_url)

links <- webpage |> 
  html_nodes("a") |>
  html_attr("href")

zip_links <- links |> 
  as_tibble() |> 
  filter(str_detect(value, "\\.zip$")) |> 
  mutate(full_url = paste0(base_url, value)) |> # Construct full URL
  pull(full_url)

download_dir <- "data/metsakeskus/metsavarakuviot"

dir.create(download_dir, showWarnings = F)

paths_in <- character()

for (url in zip_links) {
  # Extract the filename from the URL
  filename <- basename(url)
  # Define the full path for saving the file
  dest_file <- file.path(download_dir, filename)
  
  # Define the full path for saving the file
  zip_download_path <- file.path(download_dir, filename)
  paths_in[url] <- zip_download_path
  tryCatch({
    download.file(url, destfile = dest_file, mode = "wb") # 'wb' for binary files
  }, error = function(e) {
    message(paste0("Error downloading ", filename, ": ", e$message))
  })
}

