# This script extracts, summarizes, and visualizes predictions over
# environmental gradients from the fitted Joint Species Distribution Models (J-SDMs).

# INPUT:
# 1. Configuration file: "code/config_model.R" (defines MCMC parameters and runs).
# 2. Fitted Hmsc model objects (".rds" format) located in "models/[run_name]/".
# 3. Covariate dictionary: "data/covariates/mapping_covariate_functions.xlsx"

# OUTPUT:
# Outputs are generated dynamically for each model run and stored in the 
# respective directory: "models/[run_name]/gradient_simulations/".
#
# Generated files include:
# - [covariate]_gradient_simulations.pdf: A single PDF per covariate containing 
#   visualizations of summed responses, example species, and community-weighted mean traits.

# Source configuration
source("code/config_model.R")

library(Hmsc)
library(ggplot2)
library(readxl)
library(dplyr)
library(here)

setwd(here::here())

# Set default parameters for prediction evaluation
species_list <- NULL # Default: one example species shown, selected as prevalence closest to 0.5 (probit) or most abundant (other)
trait_list <- NULL   # Default: community weighted mean shown for all traits

# Iterate through each model configuration
for (i in 1:nrow(run_config$mcmc)) {
  # i <- 1
  run_name <- generate_run_name(run_config)[i]
  
  # Define paths
  model_dir <- file.path("models", run_name)
  responses_dir <- file.path(model_dir, "gradient_simulations")
  
  if (!dir.exists(responses_dir)) {
    dir.create(responses_dir, recursive = TRUE)
  }
  
  model_path <- file.path(model_dir, paste0("fitted_", run_name, ".rds"))
  
  if (file.exists(model_path)) {
    message(paste("Visualizing responses for run:", run_name))
    
    m <- readRDS(model_path)
    
    # Reconstruct the orthogonal polynomial bases from the original training data
    poly_last <- poly(m$XData$mean_temperature_yr_last, degree = 2)
    poly_yr <- poly(m$XData$mean_temperature_yr, degree = 2)
    
    # Create a temporary model object that bypasses internal formula evaluation
    m_temp <- m
    m_temp$XFormula <- as.formula("~.")
    
    # Read the covariate dictionary and filter variables for gradient simulations
    dict <- readxl::read_excel("data/covariates/mapping_covariate_functions.xlsx")
    cov_names_model <- m$covNames
    
    cov_map <- data.frame(model_cov = cov_names_model, var_target = cov_names_model) |>
      dplyr::left_join(dplyr::select(dict, var_target, type), by = "var_target")
    
    # Select only covariates marked with 1 in 'to_gradient'
    covariates <- cov_map$model_cov[!is.na(cov_map$type) & cov_map$type == "forest_structure"]
    
    # Determine the example species
    ex_sp <- which.max(colMeans(m$Y, na.rm = TRUE))
    if (m$distr[1,1] == 2) {
      ex_sp <- which.min(abs(colMeans(m$Y, na.rm = TRUE) - 0.5))
    }
    if (!is.null(species_list)) {
      ex_sp <- species_list
    }
    
    # Loop through filtered covariates to construct gradients and plot
    if (length(covariates) > 0) {
      for (k in 1:(length(covariates))) {
        # k <- 1
        covariate <- covariates[[k]]
        
        Gradient <- constructGradient(m, focalVariable = covariate)
        Gradient2 <- constructGradient(m, focalVariable = covariate, non.focalVariables = 1)
        
        # Inject the exact polynomial terms using the training basis
        inject_poly <- function(grad_data) {
          p_last <- predict(poly_last, grad_data$mean_temperature_yr_last)
          p_yr <- predict(poly_yr, grad_data$mean_temperature_yr)
          
          grad_data$`poly(mean_temperature_yr_last, degree = 2)1` <- p_last[, 1]
          grad_data$`poly(mean_temperature_yr_last, degree = 2)2` <- p_last[, 2]
          grad_data$`poly(mean_temperature_yr, degree = 2)1` <- p_yr[, 1]
          grad_data$`poly(mean_temperature_yr, degree = 2)2` <- p_yr[, 2]
          
          return(grad_data)
        }
        
        Gradient$XDataNew <- inject_poly(Gradient$XDataNew)
        Gradient2$XDataNew <- inject_poly(Gradient2$XDataNew)
        
        # Predict values. The output is structured as a three-dimensional array: 
        # sites (rows) x species (columns) x posterior samples (third dimension).
        predY <- predict(m_temp, Gradient = Gradient, expected = TRUE)  
        predY2 <- predict(m_temp, Gradient = Gradient2, expected = TRUE)  
        
        # Open PDF device for the specific covariate
        pdf(file = file.path(responses_dir, paste0(covariate, "_simul.pdf")))
        
        # 1. Summed Response Plots
        par(mfrow = c(2, 1))
        pl <- plotGradient(m, Gradient, pred = predY, yshow = 0, measure = "S", showData = TRUE, 
                           main = paste0(run_name, ": summed response (total effect)"))
        if (inherits(pl, "ggplot")) {
          print(pl + labs(title = paste0(run_name, ": summed response (total effect)")))
        }
        
        pl <- plotGradient(m, Gradient2, pred = predY2, yshow = 0, measure = "S", showData = TRUE, 
                           main = paste0(run_name, ": summed response (marginal effect)"))
        if (inherits(pl, "ggplot")) {
          print(pl + labs(title = paste0(run_name, ": summed response (marginal effect)")))
        }
        
        # 2. Example Species Plots
        for (l in 1:length(ex_sp)) {
          par(mfrow = c(2, 1))
          pl <- plotGradient(m, Gradient, pred = predY, yshow = if(m$distr[1,1] == 2){c(-0.1, 1.1)} else {0}, 
                             measure = "Y", index = ex_sp[l], showData = TRUE, 
                             main = paste0(run_name, ": example species (total effect)"))
          if (inherits(pl, "ggplot")) {
            print(pl + labs(title = paste0(run_name, ": example species (total effect)")))
          }
          
          pl <- plotGradient(m, Gradient2, pred = predY2, yshow = if(m$distr[1,1] == 2){c(-0.1, 1.1)} else {0}, 
                             measure = "Y", index = ex_sp[l], showData = TRUE, 
                             main = paste0(run_name, ": example species (marginal effect)"))
          if (inherits(pl, "ggplot")) {
            print(pl + labs(title = paste0(run_name, ": example species (marginal effect)")))
          }
        }
        
        # 3. Community Weighted Mean Trait Plots
        if (m$nt > 1) {
          traitSelection <- 2:m$nt
          if (!is.null(trait_list)) traitSelection <- trait_list
          
          for (l in traitSelection) {
            par(mfrow = c(2, 1))
            pl <- plotGradient(m, Gradient, pred = predY, measure = "T", index = l, showData = TRUE, yshow = 0,
                               main = paste0(run_name, ": community weighted mean trait (total effect)"))
            if (inherits(pl, "ggplot")) {
              print(pl + labs(title = paste0(run_name, ": community weighted mean trait (total effect)")))
            }
            
            pl <- plotGradient(m, Gradient2, pred = predY2, measure = "T", index = l, showData = TRUE, yshow = 0,
                               main = paste0(run_name, ": community weighted mean trait (marginal effect)"))
            if (inherits(pl, "ggplot")) {
              print(pl + labs(title = paste0(run_name, ": community weighted mean trait (marginal effect)")))
            }
          }
        }
        
        # Close PDF device for the specific covariate
        dev.off()
      }
    }
    
  } else {
    warning(paste("Model file not found for run:", run_name, "- Skipping."))
  }
}