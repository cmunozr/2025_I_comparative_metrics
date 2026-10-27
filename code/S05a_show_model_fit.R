#' DESCRIPTION
#' S05a_show_model_fit.R: Iterates through the MCMC configuration grid,
#' compares full-model explanatory power metrics against cross-validation and
#' holdout predictive metrics, and generates scatter plots.
#'
#' METHODOLOGY
#' - Loads model configuration settings and utility routines.
#' - Reads explanatory model fit (mf_*.rds) and evaluation fit (mfeval_*.rds) files.
#' - Dynamically detects model types and evaluates performance metrics:
#'   - Presence-Absence: Tjur R2 and AUC.
#'   - Continuous / Count Abundance: Spearman's R2 (SR2).
#'   - Both types: Root Mean Squared Error (RMSE).
#' - Evaluates validation strategies (route-blocked CV, METSO holdout, North-South holdout).
#' - Plots explanatory vs. predictive scatter plots with 1:1 reference lines.
#'
#' INPUTS
#' - code/config_model.R: Central configuration parameters (run_config, mcmc grid, model_id).
#' - code/_utilities_hmsc_gpu.R: GPU helper functions (generate_run_name).
#' - models/{run_name}/model_fit/mf_{run_name}.rds: Baseline full-model explanatory metrics.
#' - models/{run_name}/model_fit/mfeval_{run_name}_{label}_.rds: Cross-validation / holdout predictive metrics.
#'
#' OUTPUTS
#' - models/model_fit_{model_id}_validation_comparison.pdf: Multi-page scatter plot PDF comparing explanatory vs. predictive fit.

library(Hmsc)
library(jsonify)
library(here)
library(dplyr)
source(file.path("code", "config_model.R"))
source(file.path("code", "_utilities_hmsc_gpu.R"))
set.seed(110724)

# Define vector of strategies to plot
validation_strategies <- c("route_blocked_cv", "metso_holdout", "north_south")

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
    } else if (strategy == "north_south") {
      label <- "north_south"
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
