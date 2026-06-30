#' @title Model Fit Diagnostics
#'
#' @description
#' Iterates through the compiled MCMC configuration grid and extracts both 
#' full-model explanatory power metrics and cross-validation/hold-out predictive 
#' metrics. For each model version and validation strategy, it generates 
#' performance diagnostic scatter plots mapping explanatory power against predictive 
#' power across all modeled species.
#'
#' @details
#' The script reads baseline explanatory metrics (\code{mf_*.rds}) and cross-validation 
#' evaluation metrics (\code{mfeval_*.rds}) from the designated model fit directory. 
#' It dynamically adjusts the plotting metrics based on the data structure detected in 
#' the saved objects:
#' \itemize{
#'   \item \strong{Presence-Absence Models (PA):} Quantified via Tjur $R^2$ and AUC.
#'   \item \strong{Continuous Abundance Models:} Quantified via Spearman's $R^2$ (SR2).
#'   \item RMSE is calculated for both kind of model (presence-abscence or continuous)
#' }
#' 
#' @section Validation Strategies Covered:
#' \itemize{
#'   \item \code{route_blocked_cv}: Spatial cross-validation blocking at the route level to handle spatial autocorrelation.
#'   \item \code{metso_holdout}: Non-random domain hold-out evaluating extrapolation to conservation-selected contexts.
#'   \item \code{random_cv}: Completely randomized sample-unit cross-validation.
#' }
#'
#'
#' @return The script writes a compiled PDF file named \code{model_fit_[model_id].pdf} 
#' into the primary \code{models} directory.
#'

library(Hmsc)
library(jsonify)
library(here)
library(dplyr)
source(file.path("code", "config_model.R"))
source(file.path("code", "_utilities_hmsc_gpu.R"))
set.seed(110724)

# Define vector of strategies to plot
validation_strategies <- c("route_blocked_cv", "metso_holdout", "random_cv")

# --- 2. Configuration and Setup ---
models_dir <- file.path(here::here(), "models")

pdf_output_path <- file.path(models_dir, paste0("model_fit_", run_config$model_id, "_validation_comparison.pdf"))
pdf(file = pdf_output_path, width = 8, height = 8)

message(paste0("--- Generating Performance Plots for Model ID: ", run_config$model_id))

# --- 3. Main Loop: Iterate Over MCMC Configurations ---
for (i in 1:nrow(run_config$mcmc)) {
  run_name <- generate_run_name(run_config)[i]
  base_model_name <- run_name 
  dir_fit <- file.path(models_dir, base_model_name, "model_fit")
  
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  
  # Load the explanatory performance file (Full Model)
  expl_output_path <- file.path(dir_fit, paste0("mf_", base_model_name, ".rds"))
  if (!file.exists(expl_output_path)) {
    message("  Explanatory model fit file not found. Skipping config.")
    next
  }
  cMF <- readRDS(expl_output_path)
  
  for (strategy in validation_strategies) {
    
    if (strategy == "metso_holdout") {
      label <- "ho_metso"
    } else if (strategy == "route_blocked_cv") {
      label <- "cv_route"
    } else if (strategy == "random_cv") {
      label <- "cv_random"
    }
    
    mf_output_path <- file.path(dir_fit, paste0("mfeval_", base_model_name, "_", label, "_", ".rds"))
    
    if (!file.exists(mf_output_path)) {
      message("    Evaluation file missing for ", label, ". Skipping strategy.")
      next
    }
    cMFCV <- readRDS(mf_output_path)
    
    message("    Plotting strategy: ", strategy)
    
    # --- Plotting ---
    
    # Tjur R2 (Presence-Absence)
    if (!is.null(cMF$TjurR2) && !is.null(cMFCV$TjurR2)) {
      plot(cMF$TjurR2, cMFCV$TjurR2, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "Explanatory Power",
           ylab = "Predictive Power",
           main = paste0(run_name, "\nStrategy: ", strategy, " | Tjur R2\n",
                         "mean(MF) = ", round(mean(cMF$TjurR2, na.rm = TRUE), 3),
                         ", mean(MFCV) = ", round(mean(cMFCV$TjurR2, na.rm = TRUE), 3)))
      abline(0, 1, col = "red", lty = 2)
      abline(v = 0, h = 0, col = "gray")
    }
    
    # 2. RMSE
    if (!is.null(cMF$RMSE) && !is.null(cMFCV$RMSE)) {
      plot(cMF$RMSE, cMFCV$RMSE, 
           xlim = c(0, max(c(cMF$RMSE, cMFCV$RMSE), na.rm = TRUE)), 
           ylim = c(0, max(c(cMF$RMSE, cMFCV$RMSE), na.rm = TRUE)),
           xlab = "Explanatory RMSE (Full Model)",
           ylab = "Predictive RMSE (Cross-Validation)",
           main = paste0(run_name, "\nStrategy: ", strategy, " | RMSE by Species\n",
                         "mean(MF) = ", round(mean(cMF$RMSE, na.rm = TRUE), 3),
                         ", mean(MFCV) = ", round(mean(cMFCV$RMSE, na.rm = TRUE), 3)))
      abline(0, 1, col = "red", lty = 2)
    }
    
    # AUC (Presence-Absence discriminating capability)
    if (!is.null(cMF$AUC) && !is.null(cMFCV$AUC)) {
      plot(cMF$AUC, cMFCV$AUC, xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Explanatory Power",
           ylab = "Predictive Power",
           main = paste0(run_name, "\nStrategy: ", strategy, " | AUC\n",
                         "mean(MF) = ", round(mean(cMF$AUC, na.rm = TRUE), 3),
                         ", mean(MFCV) = ", round(mean(cMFCV$AUC, na.rm = TRUE), 3)))
      abline(0, 1, col = "red", lty = 2)
      abline(v = 0.5, h = 0.5, col = "gray")
    }
    
    # Spearmans R2 (SR2 - Continuous Abundance / Counts)
    if (!is.null(cMF$SR2) && !is.null(cMFCV$SR2)) {
      plot(cMF$SR2, cMFCV$SR2, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "Explanatory Power",
           ylab = "Predictive Power",
           main = paste0(run_name, "\nStrategy: ", strategy, " | SR2\n",
                         "mean(MF) = ", round(mean(cMF$SR2, na.rm = TRUE), 3),
                         ", mean(MFCV) = ", round(mean(cMFCV$SR2, na.rm = TRUE), 3)))
      abline(0, 1, col = "red", lty = 2)
      abline(v = 0, h = 0, col = "gray")
    }   
  }
}

dev.off()
