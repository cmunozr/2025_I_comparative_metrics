# --- 1. Load Libraries ---
library(tidyverse) # Includes dplyr, stringr, purrr, ggplot2
library(colorspace)

# --- 2. Configuration ---
# Define where your models live
base_dir <- "models"

# Ensure output directory exists
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. This function takes a folder path, reads ALL files, and returns a list
process_model_folder <- function(folder_path) {
  
  run_name <- basename(folder_path)
  fit_dir  <- file.path(folder_path, "model_fit")
  
  # 3.1. Validation
  if (!dir.exists(fit_dir)) {
    message(paste("Skipping (no fit folder):", run_name))
    return(NULL)
  }
  
  message(paste("Processing:", run_name))
  
  # 3.2. Locate Files
  files <- list.files(fit_dir, full.names = TRUE)
  
  path_mf     <- files[str_detect(files, "mf_")][1]
  path_mfeval <- files[str_detect(files, "mfeval")][1]
  path_waic   <- files[str_detect(files, "waic")][1]
  

  # 3.3. Extract Metadata from Folder Name (Regex)
  # We extract these now so they are available in both dataframes
  
  short_model_names <- stringr::str_extract(run_name, "fbs_M[0-9]{3}(\\.[1-9])?")
  
  meta_thin    <- str_extract(run_name, "thin_[0-9]+")
  meta_samples <- str_extract(run_name, "samples_[0-9]+")
  meta_chains  <- str_extract(run_name, "chains_[0-9]+")
  
  # Create a unique class label
  model_class  <- paste(meta_thin, meta_samples, meta_chains, sep = "_")
  
  # 3.4. Read Data
  waic_val   <- read_rds(path_waic) %>% unlist() %>% as.numeric()
  mf_raw     <- read_rds(path_mf)
  mfeval_raw <- read_rds(path_mfeval)
  
  # 3.5. Build mf_df part
  # Combine list elements, convert to DF, add columns
  df_mf <- do.call("cbind", mf_raw) %>% 
    as.data.frame() %>%
    mutate(
      Model = short_model_names,
      WAIC = waic_val,
      Class = model_class,
      Type = "Fit_Stats"
    )
  
  # 3.6. Build mfeval_df part
  df_eval <- do.call("cbind", mfeval_raw) %>% 
    as.data.frame() %>%
    mutate(
      Model = short_model_names,
      WAIC = waic_val, # WAIC is often useful here too
      Class = model_class,
      Type = "Eval_Stats"
    )
  
  # Return both as a named list
  return(list(mf = df_mf, mfeval = df_eval))
}

# --- 4. Execution ---

message("--- Starting Data Aggregation ---")

# Get list of directories
model_dirs <- list.dirs(base_dir, recursive = FALSE) %>% 
  str_subset("fbs_M")

# Run the function using purrr::map
# This creates a large list of lists
all_results <- map(model_dirs, process_model_folder)

# Remove any NULL results (failed folders)
all_results <- compact(all_results)

# --- 5. Split and Bind ---

# Extract all 'mf' items and bind them
mf_df <- map(all_results, "mf") %>% 
  bind_rows()

# Extract all 'mfeval' items and bind them
mfeval_df <- map(all_results, "mfeval") %>% 
  bind_rows()

# reorder WAIC

mf_df <- mf_df %>%
  mutate(WAIC_Label = fct_reorder(as.factor(round(WAIC, 1)), WAIC))
mfeval_df <- mfeval_df %>%
  mutate(WAIC_Label = fct_reorder(as.factor(round(WAIC, 1)), WAIC))

# --- 6. Plotting A---

types <- c("mf_df", "mfeval_df")
stats <- c("SR2", "RMSE")
lims <- list(c(-0.25,1), c(0,4))


for(t in 1:length(types)){
  types_t <- types[t]
  for(s in 1:length(stats)){
    stats_s <- stats[s]
    model_comparison_plot <- ggplot(get(types_t), aes(x = WAIC_Label, y = get(stats_s), fill = Model)) +
      geom_boxplot(width = 0.3, alpha = 0.5) +
      stat_summary(
        fun = mean, geom = "text", 
        aes(label = paste0("\u03BC=", round(after_stat(y), 3))),
        vjust = -0.8,   # Moves text slightly up
        size = 3,       # Text size
        color = "black", fontface = "bold"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      coord_cartesian(ylim = lims[[s]], clip = "off") +
      labs(
        title = paste0(stats_s, " vs. Accuracy (WAIC)"),
        subtitle = "Model Comparison",
        x = "WAIC (Lower is better)",
        y = stats_s
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    f <- file.path("models", paste0(types_t, "_", stats_s, ".jpeg"))
    ggplot2::ggsave(filename = f, model_comparison_plot, width = 7, height = 7)
  }
}




