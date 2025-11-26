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
convergence_folder[sapply(convergence_folder, file.exists)]

output_dir <- file.path("models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message(paste0("--- Starting comparison for all models runned ---"))

# Initialize lists to store aggregated data from all runs
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
for (i in 1:length(convergence_folder)) {
  
  if(is.null(sum_func)){
    stop("Provide an object function to summarize")
  }else{
    if(is.null(sum_name)){
      stop("Provide a name to the function")
    }
  }
  
  run_name <- model_setup[i]
  evaluation_dir <- convergence_folder[i]
  
  if (!dir.exists(evaluation_dir)) {
    message(paste("Evaluation directory for", run_name, "not found. Skipping."))
    next
  }
  
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

run_names_found <- names(psrf.beta.all) 

short_model_names <- stringr::str_extract(run_names_found, "fbs_M[0-9]{3}(\\.[1-9])?")

thin_vals <- stringr::str_extract(run_names_found, "thin_[0-9]+")
sample_vals <- stringr::str_extract(run_names_found, "samples_[0-9]+")
chain_vals <- stringr::str_extract(run_names_found, "chains_[0-9]+")
model_classes <- paste(thin_vals, sample_vals, chain_vals, sep = "_")
unique_classes <- unique(model_classes)
class_shape_legend_map <- setNames(1:length(unique_classes), unique_classes)
model_shapes_for_plot <- class_shape_legend_map[model_classes]
names(model_shapes_for_plot) <- run_names_found

ess_min_threshold <- 1000

calc_stat <- function(data_list, stat_func = mean) {
  sapply(data_list, stat_func, na.rm = TRUE)
}
ess_param_lists <- list(
  "Beta" = ess.beta.all,
  "Gamma" = ess.gamma.all,
  "Rho" = ess.rho.all,
  "V" = ess.V.all,
  "Omega_vakio" = ess.omega.vakio.all,
  "Omega_year" = ess.omega.year.all,
  "Omega_sampleUnit" = ess.omega.sampleUnit.all,
  "Alpha" = ess.alpha.all
)
psrf_param_lists <- list(
  "Beta" = psrf.beta.all,
  "Gamma" = psrf.gamma.all,
  "Rho" = psrf.rho.all,
  "V" = psrf.V.all,
  "Omega_vakio" = psrf.omega.vakio.all,
  "Omega_year" = psrf.omega.year.all,
  "Omega_sampleUnit" = psrf.omega.sampleUnit.all,
  "Alpha" = psrf.alpha.all
)

mean_ess_summary <- data.frame(lapply(ess_param_lists, calc_stat, stat_func = sum_func))
mean_psrf_summary <- data.frame(lapply(psrf_param_lists, calc_stat, stat_func = sum_func))

rownames(mean_ess_summary) <- run_names_found
rownames(mean_psrf_summary) <- run_names_found

write.csv(mean_ess_summary, 
          file = file.path(output_dir, "fbs_Models_all_ess_summary.csv"), 
          row.names = TRUE)
write.csv(mean_psrf_summary, 
          file = file.path(output_dir, "fbs_Models_all_psrf_summary.csv"), 
          row.names = TRUE)


# --- 5. Plotting Section ---

plot_convergence_dot_grid <- function(ess_data_list, psrf_data_list, param_name, 
                                      shape_vector, legend_map, short_names,
                                      stat_func = mean, stat_name = "Mean") { 
  
  if(length(ess_data_list) == 0 || length(psrf_data_list) == 0) {
    message(paste("Skipping", param_name, "- missing ESS or PSRF data."))
    return()
  }
  
  ess_values <- sapply(ess_data_list, stat_func, na.rm = TRUE)
  psrf_values <- sapply(psrf_data_list, stat_func, na.rm = TRUE)
  
  model_names <- names(shape_vector) 
  ess_values <- ess_values[model_names]
  psrf_values <- psrf_values[model_names]
  
  short_names_ordered <- short_names[match(model_names, run_names_found)]
  x_pos <- 1:length(model_names)
  
  ess_max_full <- max(c(0, unlist(ess_data_list)), na.rm = TRUE)
  ess_ylim_full <- c(0, ess_max_full)
  ess_ylim_zoom <- c(0, ess_min_threshold * 3) 
  
  psrf_max_full <- max(c(1.2, unlist(psrf_data_list)), na.rm = TRUE)
  psrf_ylim_full <- c(0.9, psrf_max_full)
  psrf_ylim_zoom <- c(0.98, 1.1)
  
  par(mfrow = c(2, 2))
  par(oma = c(0, 0, 3, 0)) 
  
  # --- PLOT 1: ESS (Full) ---
  par(mar = c(3, 4, 4, 2) + 0.1) 
  plot(x_pos, ess_values,
       ylim = ess_ylim_full,
       main = "ESS (Full)",
       ylab = paste(stat_name, "ESS"),
       xlab = "",
       xaxt = "n", 
       col = "black",
       pch = shape_vector[model_names],
       cex = 1.5)
  abline(h = ess_min_threshold, col = "red", lty = 2)
  legend("topright", legend = names(legend_map), pch = legend_map, 
         col = "black", cex = 0.6, bty = "n")
  
  # --- PLOT 2: ESS (Zoomed) ---
  par(mar = c(3, 4, 4, 2) + 0.1)
  plot(x_pos, ess_values,
       ylim = ess_ylim_zoom,
       main = "ESS (Zoomed)",
       ylab = paste(stat_name, "ESS"),
       xlab = "",
       xaxt = "n",
       col = "black",
       pch = shape_vector[model_names],
       cex = 1.5)
  abline(h = ess_min_threshold, col = "red", lty = 2)
  
  # --- PLOT 3: PSRF (Full) ---
  par(mar = c(8, 4, 4, 2) + 0.1) 
  plot(x_pos, psrf_values,
       ylim = psrf_ylim_full,
       main = "PSRF (Full)",
       ylab = paste(stat_name, "PSRF"),
       xlab = "",
       xaxt = "n",
       col = "black",
       pch = shape_vector[model_names],
       cex = 1.5)
  abline(h = 1.1, col = "red", lty = 2)
  abline(h = 1.05, col = "green", lty = 2)
  axis(1, at = x_pos, labels = short_names_ordered, las = 2, cex.axis = 0.8)
  
  # --- PLOT 4: PSRF (Zoomed) ---
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot(x_pos, psrf_values,
       ylim = psrf_ylim_zoom,
       main = "PSRF (Zoomed)",
       ylab = paste(stat_name, "PSRF"),
       xlab = "",
       xaxt = "n",
       col = "black",
       pch = shape_vector[model_names],
       cex = 1.5)
  abline(h = 1.1, col = "red", lty = 2)
  abline(h = 1.05, col = "green", lty = 2)
  axis(1, at = x_pos, labels = short_names_ordered, las = 2, cex.axis = 0.8)
  
  mtext(paste("Parameter:", param_name, "(", stat_name, ")"), 
        outer = TRUE, 
        cex = 1.5, 
        font = 2, 
        line = 1)
}


# --- PDF GENERATION ---
pdf(file = file.path(output_dir, "fbs_Models_all_psrf_ess.pdf"), width = 11, height = 8.5)

plot_convergence_dot_grid(ess.beta.all, psrf.beta.all, "Beta", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.gamma.all, psrf.gamma.all, "Gamma", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.rho.all, psrf.rho.all, "Rho", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.V.all, psrf.V.all, "V", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.omega.vakio.all, psrf.omega.vakio.all, "Omega (vakio)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.omega.year.all, psrf.omega.year.all, "Omega (year)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.omega.sampleUnit.all, psrf.omega.sampleUnit.all, "Omega (sampleUnit)", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)
plot_convergence_dot_grid(ess.alpha.all, psrf.alpha.all, "Alpha", 
                          model_shapes_for_plot, class_shape_legend_map, short_model_names, sum_func, sum_name)

dev.off()
