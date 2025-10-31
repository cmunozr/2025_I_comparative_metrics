# This script automatically evaluates the MCMC convergence for EACH run defined
# in the config_model.R grid. For each run, it creates a dedicated 'convergence'
# sub-folder containing:
#   1. A text report ('convergence_report.txt').
#   2. A PDF of MCMC trace plots for key parameters ('trace_plots.pdf').
#   3. Raw data CSVs for PSRF and ESS values, to be used by the S03b.

# --- 1. Load Libraries and Configuration ---
library(Hmsc)
library(colorspace)
library(coda)
library(tidyverse)

set.seed(11072024)
source("code/config_model.R")

# --- 2. Configuration for which parameters to evaluate ---
showBeta <- TRUE
showGamma <- TRUE
showOmega <- TRUE
maxOmega <- 100
showRho <- TRUE
showAlpha <- TRUE
showV <- TRUE
num_trace_plots <- 5

# --- 3. Main Loop: Iterate Over the MCMC Grid ---
for (i in 1:nrow(run_config$mcmc)) {
  # i <- 3
  # --- A. Setup for the Current Run ---
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  
  fitted_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  evaluation_dir <- file.path("models", run_name, "convergence")
  
  # --- B. Pre-run Checks ---
  if (!file.exists(fitted_model_path)) {
    message("Fitted model file not found. Skipping.")
    next
  }
  
  report_file <- file.path(evaluation_dir, "convergence_report.txt")
  if (file.exists(report_file)) {
    message("Convergence report already exists. Skipping.")
    next
  }
  
  # --- C. Evaluation of convergence ---
  tryCatch({
    dir.create(evaluation_dir, recursive = TRUE, showWarnings = FALSE)
    fitted_model <- readRDS(fitted_model_path)
    
    if(run_config$mcmc$n_chains[i] < 2) {
      warning("Convergence diagnostics require at least 2 chains. Skipping run.")
      next
    }
    
    # Open PDF for this run's trace plots
    trace_plots_file <- file.path(evaluation_dir, "trace_plots.pdf")
    pdf(file = trace_plots_file)
    
    # Start the text report
    cat(paste0("MCMC Convergence Statistics for run: ", run_name, "\n\n"), file = report_file)
    
    mpost <- convertToCodaObject(fitted_model, spNamesNumbers = c(TRUE,FALSE), covNamesNumbers = c(TRUE,FALSE))
    sum_text <- "Min.,  1st Qu.,  Median,  Mean,  3rd Qu.,  Max"
    
    # --- BETA ---
    if(showBeta){
      cat("\n--- Beta (Species Niches) ---\n", file=report_file, sep="", append=TRUE)
      psrf <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
      ess <- effectiveSize(mpost$Beta)
      
      cat("\nPSRF Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(psrf[,1]), file=report_file, append=TRUE)
      cat("\nESS Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(ess), file=report_file, append=TRUE)
      
      write.csv(data.frame(psrf), file.path(evaluation_dir, "psrf_beta.csv"))
      write.csv(data.frame(ess), file.path(evaluation_dir, "ess_beta.csv"))
      
      # Trace Plots
      par_names <- colnames(mpost$Beta[[1]]); num_to_plot <- min(num_trace_plots, length(par_names))
      if(num_to_plot > 0) {
        selected_indices <- sample(1:length(par_names), num_to_plot)
        for(p_name in par_names[selected_indices]) { plot(mpost$Beta[, p_name, drop=FALSE], main=paste("Trace of Beta:", p_name)) }
      }
    }
    
    # --- GAMMA ---
    if(showGamma){
      cat("\n\n--- Gamma (Traits) ---\n", file=report_file, sep="", append=TRUE)
      psrf <- gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
      ess <- effectiveSize(mpost$Gamma)
      
      cat("\nPSRF Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(psrf[,1]), file=report_file, append=TRUE)
      cat("\nESS Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(ess), file=report_file, append=TRUE)
      
      write.csv(data.frame(psrf), file.path(evaluation_dir, "psrf_gamma.csv"))
      write.csv(data.frame(ess), file.path(evaluation_dir, "ess_gamma.csv"))
      
      # Trace Plots
      par_names <- colnames(mpost$Gamma[[1]]); num_to_plot <- min(num_trace_plots, length(par_names))
      if(num_to_plot > 0) {
        selected_indices <- sample(1:length(par_names), num_to_plot)
        for(p_name in par_names[selected_indices]) { plot(mpost$Gamma[, p_name, drop=FALSE], main=paste("Trace of Gamma:", p_name)) }
      }
    }
    
    # --- RHO (Phylogenetic Signal) ---
    if(showRho && !is.null(mpost$Rho)){
      cat("\n\n--- Rho (Phylogenetic Signal) ---\n", file=report_file, sep="", append=TRUE)
      rho <- lapply(mpost$Rho, as.matrix)
      class(rho) <- "mcmc.list"
      psrf <- coda::gelman.diag(rho, multivariate=FALSE, autoburnin = FALSE)$psrf
      ess <- effectiveSize(mpost$Rho)
      
      cat("\nPSRF:", psrf[1], file=report_file, append=TRUE)
      cat("\nESS:", ess[1], file=report_file, append=TRUE)
      
      write.csv(data.frame(psrf), file.path(evaluation_dir, "psrf_rho.csv"))
      write.csv(data.frame(ess), file.path(evaluation_dir, "ess_rho.csv"))
      
      plot(mpost$Rho, main="Trace of Rho")
    }
    
    # --- V (Covariates) ---
    if(showV && !is.null(mpost$V)){
      cat("\n\n--- V (Covariates) ---\n", file=report_file, sep="", append=TRUE)
      psrf <- gelman.diag(mpost$V, multivariate=FALSE)$psrf
      ess <- effectiveSize(mpost$V)
      
      cat("\nPSRF Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(psrf[,1]), file=report_file, append=TRUE)
      cat("\nESS Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(ess), file=report_file, append=TRUE)
      
      write.csv(data.frame(psrf), file.path(evaluation_dir, "psrf_v.csv"))
      write.csv(data.frame(ess), file.path(evaluation_dir, "ess_v.csv"))
    }
    
    # --- OMEGA (Species Associations) ---
    if(showOmega && length(mpost$Omega) > 0){
      cat("\n\n--- Omega (Species Associations) ---\n", file=report_file, sep="", append=TRUE)
      for(k in 1:length(mpost$Omega)){
        rl_name <- names(fitted_model$ranLevels)[k]
        cat(paste0("\nRandom Level: ", rl_name, "\n"), file=report_file, append=TRUE)
        
        tmp <- mpost$Omega[[k]]
        # Sub-sample if matrix is very large
        if(ncol(tmp[[1]]) > maxOmega^2){
          selected_cols <- sample(1:ncol(tmp[[1]]), maxOmega^2)
          tmp <- tmp[, selected_cols, drop=FALSE]
        }
        
        psrf <- gelman.diag(tmp, multivariate=FALSE)$psrf
        ess <- effectiveSize(tmp)
        
        cat("\nPSRF Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(psrf[,1]), file=report_file, append=TRUE)
        cat("\nESS Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(ess), file=report_file, append=TRUE)
        
        write.csv(data.frame(psrf), file.path(evaluation_dir, paste0("psrf_omega_", rl_name, ".csv")))
        write.csv(data.frame(ess), file.path(evaluation_dir, paste0("ess_omega_", rl_name, ".csv")))
        
        # Trace Plots
        par_names <- colnames(tmp[[1]]); num_to_plot <- min(num_trace_plots, length(par_names))
        if(num_to_plot > 0) {
          selected_indices <- sample(1:length(par_names), num_to_plot)
          for(p_name in par_names[selected_indices]) { plot(tmp[, p_name, drop=FALSE], main=paste("Trace of Omega", rl_name, ":", p_name)) }
        }
      }
    }
    
    # --- ALPHA (Spatial Scale) ---
    if(showAlpha && length(mpost$Alpha) > 0){
      cat("\n\n--- Alpha (Spatial Autocorrelation) ---\n", file=report_file, sep="", append=TRUE)
      for(k in 1:length(mpost$Alpha)){
        rl_name <- names(fitted_model$ranLevels)[k]
        if(fitted_model$ranLevels[[k]]$sDim > 0) { # Only for spatial random levels
          cat(paste0("\nRandom Level: ", rl_name, "\n"), file=report_file, append=TRUE)
          psrf <- gelman.diag(mpost$Alpha[[k]], multivariate=FALSE)$psrf
          ess <- effectiveSize(mpost$Alpha[[k]])
          
          cat("\nPSRF Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(psrf[,1]), file=report_file, append=TRUE)
          cat("\nESS Summary:\n", file=report_file, append=TRUE); cat(sum_text, "\n", file=report_file, append=TRUE); capture.output(summary(ess), file=report_file, append=TRUE)
          
          write.csv(data.frame(psrf), file.path(evaluation_dir, paste0("psrf_alpha_", rl_name, ".csv")))
          write.csv(data.frame(ess), file.path(evaluation_dir, paste0("ess_alpha_", rl_name, ".csv")))
          
          # Trace Plots
          par_names <- colnames(mpost$Alpha[[k]][[1]]); num_to_plot <- min(num_trace_plots, length(par_names))
          if(num_to_plot > 0) {
            selected_indices <- sample(1:length(par_names), num_to_plot)
            for(p_name in par_names[selected_indices]) { plot(mpost$Alpha[[k]][, p_name, drop=FALSE], main=paste("Trace of Alpha", rl_name, ":", p_name)) }
          }
        }
      }
    }
    
    dev.off()

  }, error = function(e) {
    warning(paste("An error occurred while evaluating", run_name, ":", e$message))
    if(dev.cur() != 1) dev.off() # Ensure PDF device is closed on error
  })
}
