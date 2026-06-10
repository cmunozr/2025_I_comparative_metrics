# --- 1. Load Libraries and Configuration ---
library(tidyverse)
library(dplyr)
library(stringr)
library(colorspace) 

#setwd("E:/2025_I_comparative_metrics/")

# --- 2. Setup Directories and Initialize Data Lists ---
model_setup <- list.dirs("models", recursive = FALSE) |>
  basename() |>
  stringr::str_subset("^fbs_M")

convergence_folder <- file.path("models", model_setup, "convergence")
# Only keep paths that actually exist
existing_folders_idx <- sapply(convergence_folder, file.exists)
model_setup <- model_setup[existing_folders_idx]
convergence_folder <- convergence_folder[existing_folders_idx]

output_dir <- file.path("models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ess_summary_file <- file.path(output_dir, "fbs_Models_all_ess_summary.csv")
psrf_summary_file <- file.path(output_dir, "fbs_Models_all_psrf_summary.csv")

# Read the csvs and extract what models were evaluated
evaluated_models <- c()
if (file.exists(ess_summary_file) && file.exists(psrf_summary_file)) {
  existing_ess_summary <- read.csv(ess_summary_file, row.names = 1)
  existing_psrf_summary <- read.csv(psrf_summary_file, row.names = 1)
  evaluated_models <- rownames(existing_ess_summary)
  message(paste("Found", length(evaluated_models), "previously evaluated models in summary files."))
} else {
  existing_ess_summary <- data.frame()
  existing_psrf_summary <- data.frame()
}

# Read the folders and discard the ones that were evaluated before
models_to_evaluate_idx <- which(!model_setup %in% evaluated_models)
model_setup_new <- model_setup[models_to_evaluate_idx]
convergence_folder_new <- convergence_folder[models_to_evaluate_idx]

message(paste0("--- Processing ", length(model_setup_new), " NEW models ---"))

# Initialize lists to store aggregated data from NEW runs only
psrf.beta.all <- list(); ess.beta.all <- list()
psrf.gamma.all <- list(); ess.gamma.all <- list()
psrf.rho.all <- list(); ess.rho.all <- list()
psrf.V.all <- list(); ess.V.all <- list()
psrf.alpha.all <- list(); ess.alpha.all <- list() 
psrf.omega.vakio.all <- list(); ess.omega.vakio.all <- list()
psrf.omega.year.all <- list(); ess.omega.year.all <- list()
psrf.omega.sampleUnit.all <- list(); ess.omega.sampleUnit.all <- list()

# statistical function to summarize
sum_func <- median
sum_name <- "Median"

# --- 3. Data Aggregation Loop ---
# Process the models that are not in the summary
for (i in seq_along(convergence_folder_new)) {
  
  if(is.null(sum_func)){
    stop("Provide an object function to summarize")
  }else{
    if(is.null(sum_name)){
      stop("Provide a name to the function")
    }
  }
  
  run_name <- model_setup_new[i]
  evaluation_dir <- convergence_folder_new[i]
  
  message(paste("Aggregating data from:", run_name))
  
  read_and_store <- function(param_name, metric, data_list) {
    file_path <- file.path(evaluation_dir, paste0(metric, "_", param_name, ".csv"))
    if (file.exists(file_path)) {
      if (metric == "psrf") {
        data_list[[run_name]] <- read.csv(file_path)[, "Point.est."]
      } else { 
        data_list[[run_name]] <- read.csv(file_path)[, "ess"]
      }
    }
    return(data_list)
  }
  
  psrf.beta.all <- read_and_store("beta", "psrf", psrf.beta.all)
  ess.beta.all  <- read_and_store("beta", "ess", ess.beta.all)
  
  psrf.gamma.all <- read_and_store("gamma", "psrf", psrf.gamma.all)
  ess.gamma.all  <- read_and_store("gamma", "ess", ess.gamma.all)
  
  psrf.rho.all <- read_and_store("rho", "psrf", psrf.rho.all)
  ess.rho.all  <- read_and_store("rho", "ess", ess.rho.all)
  
  psrf.V.all <- read_and_store("v", "psrf", psrf.V.all)
  ess.V.all  <- read_and_store("v", "ess", ess.V.all)
  
  psrf.omega.vakio.all <- read_and_store("omega_vakio", "psrf", psrf.omega.vakio.all)
  ess.omega.vakio.all  <- read_and_store("omega_vakio", "ess", ess.omega.vakio.all)
  
  psrf.omega.year.all <- read_and_store("omega_year", "psrf", psrf.omega.year.all)
  ess.omega.year.all  <- read_and_store("omega_year", "ess", ess.omega.year.all)
  
  psrf.omega.sampleUnit.all <- read_and_store("omega_sampleUnit", "psrf", psrf.omega.sampleUnit.all)
  ess.omega.sampleUnit.all  <- read_and_store("omega_sampleUnit", "ess", ess.omega.sampleUnit.all)
  
  alpha_files_psrf <- list.files(evaluation_dir, pattern = "psrf_alpha_.*\\.csv", full.names = TRUE)
  if(length(alpha_files_psrf) > 0) {
    psrf.alpha.all[[run_name]] <- do.call(c, lapply(alpha_files_psrf, function(f) read.csv(f)[,"Point.est."]))
  }
  
  alpha_files_ess <- list.files(evaluation_dir, pattern = "ess_alpha_.*\\.csv", full.names = TRUE)
  if(length(alpha_files_ess) > 0) {
    ess.alpha.all[[run_name]] <- do.call(c, lapply(alpha_files_ess, function(f) read.csv(f)[,"ess"]))
  }
}

# --- 4. Calculation Section ---
run_names_found_new <- names(psrf.beta.all) 

calc_stat <- function(data_list, stat_func = mean) {
  if (length(data_list) == 0) return(NA)
  sapply(data_list, stat_func, na.rm = TRUE)
}

# Only calculate if we found new models
if (length(run_names_found_new) > 0) {
  ess_param_lists <- list(
    "Beta" = ess.beta.all, "Gamma" = ess.gamma.all, "Rho" = ess.rho.all,
    "V" = ess.V.all, "Omega_vakio" = ess.omega.vakio.all, 
    "Omega_year" = ess.omega.year.all, "Omega_sampleUnit" = ess.omega.sampleUnit.all,
    "Alpha" = ess.alpha.all
  )
  psrf_param_lists <- list(
    "Beta" = psrf.beta.all, "Gamma" = psrf.gamma.all, "Rho" = psrf.rho.all,
    "V" = psrf.V.all, "Omega_vakio" = psrf.omega.vakio.all, 
    "Omega_year" = psrf.omega.year.all, "Omega_sampleUnit" = psrf.omega.sampleUnit.all,
    "Alpha" = psrf.alpha.all
  )
  
  new_ess_summary <- data.frame(lapply(ess_param_lists, calc_stat, stat_func = sum_func))
  new_psrf_summary <- data.frame(lapply(psrf_param_lists, calc_stat, stat_func = sum_func))
  rownames(new_ess_summary) <- run_names_found_new
  rownames(new_psrf_summary) <- run_names_found_new
  
  # Append results using bind_rows to handle missing columns safely
  if (nrow(existing_ess_summary) > 0) {
    combined_ess_summary <- bind_rows(existing_ess_summary %>% rownames_to_column("model"), 
                                      new_ess_summary %>% rownames_to_column("model")) %>% 
      column_to_rownames("model")
    
    combined_psrf_summary <- bind_rows(existing_psrf_summary %>% rownames_to_column("model"), 
                                       new_psrf_summary %>% rownames_to_column("model")) %>% 
      column_to_rownames("model")
  } else {
    combined_ess_summary <- new_ess_summary
    combined_psrf_summary <- new_psrf_summary
  }
  
  # iv) Write the csvs
  write.csv(combined_ess_summary, file = ess_summary_file, row.names = TRUE)
  write.csv(combined_psrf_summary, file = psrf_summary_file, row.names = TRUE)
} else {
  message("No new models processed. Proceeding to plot existing summaries.")
  combined_ess_summary <- existing_ess_summary
  combined_psrf_summary <- existing_psrf_summary
}

# Setup plot legends based on ALL models (old + new)
all_run_names <- rownames(combined_ess_summary)
if (length(all_run_names) == 0) stop("No models available to plot.")

short_model_names <- stringr::str_remove(all_run_names, "_thin_.*$")
thin_vals <- stringr::str_extract(all_run_names, "thin_[0-9]+")
sample_vals <- stringr::str_extract(all_run_names, "samples_[0-9]+")
chain_vals <- stringr::str_extract(all_run_names, "chains_[0-9]+")
model_classes <- paste(thin_vals, sample_vals, chain_vals, sep = "_")

unique_classes <- unique(model_classes)
class_shape_legend_map <- setNames(1:length(unique_classes), unique_classes)
model_shapes_for_plot <- class_shape_legend_map[model_classes]
names(model_shapes_for_plot) <- all_run_names

ess_min_threshold <- 1000


# --- 5. Plotting Section ---

plot_convergence_dot_grid <- function(ess_values, psrf_values, param_name, 
                                      shape_vector, legend_map, short_names,
                                      stat_name = "Median") { 
  
  # Check if data exists for this parameter
  if(is.null(ess_values) || is.null(psrf_values) || all(is.na(ess_values))) {
    message(paste("Skipping", param_name, "- missing ESS or PSRF data."))
    return()
  }
  
  model_names <- names(shape_vector) 
  
  # Order the vectors based on model_names
  ess_ordered <- ess_values[model_names]
  psrf_ordered <- psrf_values[model_names]
  short_names_ordered <- short_names[match(model_names, all_run_names)]
  x_pos <- 1:length(model_names)
  
  # Use max of medians for limits, adding a 20% buffer for visual space
  ess_max_full <- max(c(0, ess_ordered), na.rm = TRUE) * 1.2
  ess_ylim_full <- c(0, ess_max_full)
  ess_ylim_zoom <- c(0, ess_min_threshold * 3) 
  
  psrf_max_full <- max(c(1.2, psrf_ordered), na.rm = TRUE)
  psrf_ylim_full <- c(0.9, psrf_max_full)
  psrf_ylim_zoom <- c(0.98, 1.1)
  
  par(mfrow = c(2, 2))
  par(oma = c(0, 0, 3, 0)) 
  
  # --- PLOT 1: ESS (Full) ---
  par(mar = c(3, 4, 4, 2) + 0.1) 
  plot(x_pos, ess_ordered, ylim = ess_ylim_full, main = "ESS (Full)",
       ylab = paste(stat_name, "ESS"), xlab = "", xaxt = "n", col = "black",
       pch = shape_vector[model_names], cex = 1.5)
  abline(h = ess_min_threshold, col = "red", lty = 2)
  legend("topright", legend = names(legend_map), pch = legend_map, 
         col = "black", cex = 0.6, bty = "n")
  
  # --- PLOT 2: ESS (Zoomed) ---
  par(mar = c(3, 4, 4, 2) + 0.1)
  plot(x_pos, ess_ordered, ylim = ess_ylim_zoom, main = "ESS (Zoomed)",
       ylab = paste(stat_name, "ESS"), xlab = "", xaxt = "n", col = "black",
       pch = shape_vector[model_names], cex = 1.5)
  abline(h = ess_min_threshold, col = "red", lty = 2)
  
  # --- PLOT 3: PSRF (Full) ---
  par(mar = c(8, 4, 4, 2) + 0.1) 
  plot(x_pos, psrf_ordered, ylim = psrf_ylim_full, main = "PSRF (Full)",
       ylab = paste(stat_name, "PSRF"), xlab = "", xaxt = "n", col = "black",
       pch = shape_vector[model_names], cex = 1.5)
  abline(h = 1.1, col = "red", lty = 2)
  abline(h = 1.05, col = "green", lty = 2)
  axis(1, at = x_pos, labels = short_names_ordered, las = 2, cex.axis = 0.8)
  
  # --- PLOT 4: PSRF (Zoomed) ---
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot(x_pos, psrf_ordered, ylim = psrf_ylim_zoom, main = "PSRF (Zoomed)",
       ylab = paste(stat_name, "PSRF"), xlab = "", xaxt = "n", col = "black",
       pch = shape_vector[model_names], cex = 1.5)
  abline(h = 1.1, col = "red", lty = 2)
  abline(h = 1.05, col = "green", lty = 2)
  axis(1, at = x_pos, labels = short_names_ordered, las = 2, cex.axis = 0.8)
  
  mtext(paste("Parameter:", param_name, "(", stat_name, ")"), 
        outer = TRUE, cex = 1.5, font = 2, line = 1)
}

# Helper to extract named vectors
get_col <- function(df, col) {
  if(col %in% colnames(df)) setNames(df[[col]], rownames(df)) else NULL
}

# --- PDF GENERATION ---
pdf_grouped_file <- file.path(output_dir, paste0("fbs_Models_all_psrf_ess.pdf"))
pdf(file = pdf_grouped_file, width = 11, height = 8.5)

plot_convergence_dot_grid(get_col(combined_ess_summary, "Beta"), get_col(combined_psrf_summary, "Beta"), "Beta", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Gamma"), get_col(combined_psrf_summary, "Gamma"), "Gamma", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Rho"), get_col(combined_psrf_summary, "Rho"), "Rho", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "V"), get_col(combined_psrf_summary, "V"), "V", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Omega_vakio"), get_col(combined_psrf_summary, "Omega_vakio"), "Omega (vakio)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Omega_year"), get_col(combined_psrf_summary, "Omega_year"), "Omega (year)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Omega_sampleUnit"), get_col(combined_psrf_summary, "Omega_sampleUnit"), "Omega (sampleUnit)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)
plot_convergence_dot_grid(get_col(combined_ess_summary, "Alpha"), get_col(combined_psrf_summary, "Alpha"), "Alpha", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_name)

dev.off()