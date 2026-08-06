# This script extracts, summarizes, and visualizes parameter estimates from the fitted 
# Joint Species Distribution Models (J-SDMs) to evaluate species responses 
# and variance partitioning.

# INPUT:
# 1. Configuration file: "code/config_model.R" (defines MCMC parameters and runs).
# 2. Fitted Hmsc model objects (".rds" format) located in "models/[run_name]/".

# OUTPUT:
# Outputs are generated dynamically for each model run and stored in the 
# respective directory: "models/[run_name]/parameter_estimates/".
#
# Generated files include:
# - parameter_estimates.pdf: Visualizations of Variance Partitioning, species 
#   niches (Beta), functional traits (Gamma), and species associations (Omega).
# - parameter_estimates_log.txt: Text log detailing parameter extraction settings.
# - parameter_estimates_VP*.csv: Numerical tables of variance partitioning results.
# - parameter_estimates_[Beta/Gamma/Omega].xlsx: Tables containing posterior means and 
#   statistical support probabilities for model parameters.

# Source configuration
source("code/config_model.R")

library(Hmsc)
library(colorspace)
library(corrplot)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(gt)

# Global switches for optional parameter evaluations
plot_gamma_switch <- FALSE
plot_omega_switch <- FALSE

# Targeted categories for Beta plot visualization
target_categories_beta <- c("Intercept", "forest_structure")

# Set default parameters for estimates evaluation
support_level_beta <- 0.95
support_level_gamma <- 0.95
support_level_omega <- 0.90
var_part_order_explained <- NULL 
var_part_order_raw <- NULL 
show_sp_names_beta <- NULL 
plot_tree <- FALSE 
omega_order <- NULL 
show_sp_names_omega <- NULL
## Uncomment to activate characteristic

#use support levels to selected the level of statistical support shown
#support_level_beta = 0.8
#support_level_gamma = 0.8
#support_level_omega = 0.8

#use var_part_order_explained to select which order species are shown in the raw variance partitioning
#var_part_order_raw should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "decreasing" if you wish to order according to explanatory power
#var_part_order_explained = list()
#var_part_order_explained[[1]] = 0
#var_part_order_explained[[2]] = c(2,1)

#use var_part_order_raw to select which order species are shown in the explained variance partitioning
#same options apply as for var_part_order_explained
#var_part_order_raw = list()
#var_part_order_raw[[1]] = "decreasing"
#var_part_order_raw[[2]] = c(1,2)

#use show_sp_names_beta to choose to show / not show species names in betaPlot
#if given, show_sp_names_beta should be a vector with length equalling number of models
#show_sp_names_beta = c(TRUE,FALSE)

#use plot_tree to choose to plot / not plot the tree in betaPlot
#if given, plot_tree should be a vector with length equalling number of models
#plot_tree = c(FALSE,FALSE)

#use omega_order to select which order species are shown in omega plots
#omega_order should be a list of length the number of models. 
#for each element provide either 0 (use default);
#or a vector of species indices;
#or "AOE" if you wish to use the angular order of the eigenvectors.
#omega_order = list()
#omega_order[[1]] = "AOE"
#omega_order[[2]] = c(2,1)
#Default: species shown in the order they are in the model
#show_sp_names_omega = c(TRUE,FALSE) #Default: species names shown in beta plot if there are at most 30 species


# Iterate through each model configuration
for (i in 1:nrow(run_config$mcmc)) {  
  run_name <- generate_run_name(run_config)[i]
  
  # Define paths
  model_dir <- file.path("models", run_name)
  estimates_dir <- file.path(model_dir, "parameter_estimates")
  
  if (!dir.exists(estimates_dir)) {
    dir.create(estimates_dir, recursive = TRUE)
  }
  
  model_path <- file.path(model_dir, paste0("fitted_", run_name, ".rds"))
  
  if (file.exists(model_path)) {
    message(paste("Extracting parameter estimates for run:", run_name))
    
    m <- readRDS(model_path)
    
    # Initialize text file for logging parameters
    text_file <- file.path(estimates_dir, "parameter_estimates_log.txt")
    cat("This file contains additional information regarding parameter estimates.\n\n", file = text_file)
    cat(c("Run: ", run_name, "\n\n"), file = text_file, sep = "", append = TRUE)
    
    pdf(file = file.path(estimates_dir, "parameter_estimates.pdf"), height = 10, width = 10)
    
    # Match model covariates to dictionary
    cov_names_model <- m$covNames
    dict <- readxl::read_excel("data/covariates/mapping_covariate_functions.xlsx")
    
    cov_map <- data.frame(model_cov = cov_names_model, var_target = cov_names_model) |>
      dplyr::left_join(dplyr::select(dict, var_target, type), by = "var_target") |>
      dplyr::mutate(type = dplyr::case_when(
        model_cov == "(Intercept)" ~ "Intercept",
        stringr::str_detect(var_target, "temp|prec") ~ "climate",
        TRUE ~ dplyr::coalesce(type, "Uncategorized")
      ))
    
    # ------------------------------------------------------------------------
    # 1. Variance Partitioning
    # ------------------------------------------------------------------------
    if (m$XFormula == "~.") {
      covariates <- colnames(m$X)[-1]
    } else {
      covariates <- attr(terms(m$XFormula), "term.labels")
    }
    
    if (m$nr + length(covariates) > 1 & m$ns > 1) {
      preds <- computePredictedValues(m)
      
      # Define groups for Hmsc variance partitioning
      group_factors <- as.factor(cov_map$type)
      group_indices <- as.numeric(group_factors)
      group_names <- levels(group_factors)
      
      VP <- computeVariancePartitioning(m, group = group_indices, groupnames = group_names)
      
      vals <- VP$vals
      mycols <- rainbow(nrow(VP$vals))
      MF <- evaluateModelFit(hM = m, predY = preds)
      R2 <- NULL
      
      if (!is.null(MF$TjurR2)) {
        TjurR2 <- MF$TjurR2
        vals <- rbind(vals, TjurR2)
        R2 <- TjurR2
      }
      if (!is.null(MF$R2)) {
        R2 <- MF$R2
        vals <- rbind(vals, R2)
      }
      if (!is.null(MF$SR2)) {
        R2 <- MF$SR2
        vals <- rbind(vals, R2)
      }
      
      filename_vp <- file.path(estimates_dir, "parameter_estimates_VP.csv")
      write.csv(vals, file = filename_vp)
      
      if (!is.null(VP$R2T$Beta)) {
        write.csv(VP$R2T$Beta, file = file.path(estimates_dir, "parameter_estimates_VP_R2T_Beta.csv"))
      }
      if (!is.null(VP$R2T$Y)) {
        write.csv(VP$R2T$Y, file = file.path(estimates_dir, "parameter_estimates_VP_R2T_Y.csv"))
      }
      
      if (is.null(var_part_order_explained) || all(var_part_order_explained == 0)) {
        c_var_part_order_explained <- 1:m$ns
      } else if (all(var_part_order_explained == "decreasing")) {
        c_var_part_order_explained <- order(R2, decreasing = TRUE)
      } else {
        c_var_part_order_explained <- var_part_order_explained
      }
      
      VP$vals <- VP$vals[, c_var_part_order_explained]
      cat(c("\nvar.part.order.explained\n\n"), file = text_file, sep = "", append = TRUE)
      cat(m$spNames[c_var_part_order_explained], file = text_file, sep = "\n", append = TRUE)
      
      plotVariancePartitioning(hM = m, VP = VP, main = paste0("Proportion of explained variance, ", run_name), 
                               cex.main = 0.8, cols = mycols, args.leg = list(bg = "white", cex = 0.7))
      
      if (is.null(var_part_order_raw) || all(var_part_order_raw == 0)) {
        c_var_part_order_raw <- 1:m$ns
      } else if (all(var_part_order_raw == "decreasing")) {
        c_var_part_order_raw <- order(R2, decreasing = TRUE)
      } else {
        c_var_part_order_raw <- var_part_order_raw
      }
      
      if (!is.null(R2)) {
        VPr <- VP
        for (k in 1:m$ns) {
          VPr$vals[, k] <- R2[k] * VPr$vals[, k]
        }
        VPr$vals <- VPr$vals[, c_var_part_order_raw]
        cat(c("\nvar.part.order.raw\n\n"), file = text_file, sep = "", append = TRUE)
        cat(m$spNames[c_var_part_order_raw], file = text_file, sep = "\n", append = TRUE)
        plotVariancePartitioning(hM = m, VP = VPr, main = paste0("Proportion of raw variance, ", run_name), 
                                 cex.main = 0.8, cols = mycols, args.leg = list(bg = "white", cex = 0.7), ylim = c(0, 1))
      }
    }
    
    # ------------------------------------------------------------------------
    # 2. Beta Parameters Evaluation & Subsetted Visualization
    # ------------------------------------------------------------------------
    if (m$nc > 1) {
      postBeta <- getPostEstimate(m, parName = "Beta")
      filename_beta <- file.path(estimates_dir, "parameter_estimates_Beta.xlsx")
      
      pivot_estimates <- function(est_matrix, value_name) {
        df <- as.data.frame(t(est_matrix))
        colnames(df) <- cov_names_model
        df$Species <- m$spNames
        tidyr::pivot_longer(df, cols = -Species, names_to = "model_cov", values_to = value_name)
      }
      
      me <- pivot_estimates(postBeta$mean, "mean")
      po <- pivot_estimates(postBeta$support, "support")
      ne <- pivot_estimates(postBeta$supportNeg, "supportNeg")
      
      beta_summary <- me |>
        dplyr::left_join(po, by = c("Species", "model_cov")) |>
        dplyr::left_join(ne, by = c("Species", "model_cov")) |>
        dplyr::left_join(cov_map, by = "model_cov") |>
        dplyr::arrange(type, model_cov, Species)
      
      vals <- list("Beta_Summary" = beta_summary)
      writexl::write_xlsx(vals, path = filename_beta)
      
      # Select targeted indices for Intercept and forest_structure
      if(!is.null(target_categories_beta)){
        selected_indices <- which(cov_map$type %in% target_categories_beta)
        selected_indices <- selected_indices[order(match(cov_map$type[selected_indices], target_categories_beta))]
      }else{
        selected_indices <- 1:length(cov_map$type)
      }
      
      
      cov_map_sub <- cov_map[selected_indices, ]
      
      postBeta_sub <- postBeta
      postBeta_sub$mean <- postBeta$mean[selected_indices, , drop = FALSE]
      postBeta_sub$support <- postBeta$support[selected_indices, , drop = FALSE]
      if (!is.null(postBeta$supportNeg)) {
        postBeta_sub$supportNeg <- postBeta$supportNeg[selected_indices, , drop = FALSE]
      }
      
      m_plot <- m
      m_plot$nc <- length(selected_indices)
      m_plot$covNames <- paste0("[", cov_map_sub$type, "] ", cov_map_sub$model_cov)
      
      c_show_sp_names <- if (is.null(show_sp_names_beta)) (is.null(m$phyloTree) && m$ns <= 30) else show_sp_names_beta
      c_plotTree <- !is.null(m$phyloTree)
      if (!is.null(plot_tree)) c_plotTree <- c_plotTree & plot_tree
      
      plotBeta(m_plot, post = postBeta_sub, supportLevel = support_level_beta, param = "Sign",
               plotTree = c_plotTree,
               covNamesNumbers = c(FALSE, TRUE),
               spNamesNumbers = c(c_show_sp_names, FALSE),
               cex = c(0.6, 0.6, 0.8))
      
      mymain <- paste0("BetaPlot (Forest Structure & Intercept), ", run_name)
      if (!is.null(m$phyloTree)) {
        mpost <- convertToCodaObject(m)
        rhovals <- unlist(poolMcmcChains(mpost$Rho))
        mymain <- paste0(mymain, ", E[rho] = ", round(mean(rhovals), 2), ", Pr[rho>0] = ", round(mean(rhovals > 0), 2))
      }
      title(main = mymain, line = 2.5, cex.main = 0.8)
      

      # reference key
      covariate_index_key <- data.frame(
        Plot_Index = paste0("C", 1:nrow(cov_map_sub)),
        Original_Covariate = cov_map_sub$model_cov,
        Category = cov_map_sub$type
      )
      
      gt_key <- covariate_index_key |>
        gt() |>
        tab_header(
          title = md("**Covariate Reference Key**"),
        ) |>
        cols_label(
          Plot_Index = "Index", Original_Covariate = "Model Covariate", Category = "Domain"
        ) |>
        cols_align(
          align = "center", columns = c(Plot_Index, Category)
        ) |>
        cols_align(
          align = "left", columns = Original_Covariate
        ) |>
        opt_row_striping() |>
        tab_options(
          table.font.size = px(12),
          heading.title.font.size = px(14),
          heading.subtitle.font.size = px(11),
          column_labels.font.weight = "bold",
          table.width = pct(100)
        )
      
      gt::gtsave(
        data = gt_key,
        filename = file.path(estimates_dir, "covariate_key_table.html")
      )
    }
    
    # ------------------------------------------------------------------------
    # 3. Gamma Parameters Evaluation (Traits)
    # ------------------------------------------------------------------------
    if (plot_gamma_switch) {
      if (!is.null(m$nt) && m$nt > 1 && m$nc > 1) {
        message("Processing Gamma parameters...")
        postGamma <- getPostEstimate(m, parName = "Gamma")
        
        plotGamma(m, post = postGamma, supportLevel = support_level_gamma, param = "Sign",
                  covNamesNumbers = c(TRUE, FALSE),
                  trNamesNumbers = c(m$nt < 21, FALSE),
                  cex = c(0.6, 0.6, 0.8))
        title(main = paste0("GammaPlot ", run_name), line = 2.5, cex.main = 0.8)
      }
    }
    
    # ------------------------------------------------------------------------
    # 4. Omega Parameters Evaluation (Species Associations)
    # ------------------------------------------------------------------------
    if (plot_omega_switch) {
      if (m$nr > 0 && m$ns > 1) {
        message("Processing Omega parameters...")
        OmegaCor <- computeAssociations(m)
        for (r in 1:m$nr) {
          toPlot <- ((OmegaCor[[r]]$support > support_level_omega) + (OmegaCor[[r]]$support < (1 - support_level_omega)) > 0) * sign(OmegaCor[[r]]$mean)
          c_show_sp_names <- if (is.null(show_sp_names_omega)) (m$ns <= 30) else show_sp_names_omega
          
          if (!c_show_sp_names) {
            colnames(toPlot) <- rep("", m$ns)
            rownames(toPlot) <- rep("", m$ns)
          }
          
          if (is.null(omega_order) || all(omega_order == 0)) {
            plotOrder <- 1:m$ns
          } else if (all(omega_order == "AOE")) {
            plotOrder <- corrMatOrder(OmegaCor[[r]]$mean, order = "AOE")
          } else {
            plotOrder <- omega_order
          }
          
          cat(c("\nomega.order\n\n"), file = text_file, sep = "", append = TRUE)
          cat(m$spNames[plotOrder], file = text_file, sep = "\n", append = TRUE)
          
          mymain <- paste0("Associations, ", run_name, ": ", names(m$ranLevels)[[r]])
          if (m$ranLevels[[r]]$sDim > 0) {
            mpost <- convertToCodaObject(m)
            alphavals <- unlist(poolMcmcChains(mpost$Alpha[[1]][, 1]))
            mymain <- paste0(mymain, ", E[alpha1] = ", round(mean(alphavals), 2), ", Pr[alpha1>0] = ", round(mean(alphavals > 0), 2))
          }
          
          corrplot(toPlot[plotOrder, plotOrder], method = "color",
                   col = colorRampPalette(c("blue", "white", "red"))(3),
                   mar = c(0, 0, 1, 0),
                   main = mymain, cex.main = 0.8)
          
          me <- as.data.frame(OmegaCor[[r]]$mean)
          me <- cbind(m$spNames, me)
          colnames(me)[1] <- ""
          po <- as.data.frame(OmegaCor[[r]]$support)
          po <- cbind(m$spNames, po)
          colnames(po)[1] <- ""
          ne <- as.data.frame(1 - OmegaCor[[r]]$support)
          ne <- cbind(m$spNames, ne)
          colnames(ne)[1] <- ""
          
          vals <- list("Posterior mean" = me, "Pr(x>0)" = po, "Pr(x<0)" = ne)
          filename_omega <- file.path(estimates_dir, paste0("parameter_estimates_Omega_", names(m$ranLevels)[[r]], ".xlsx"))
          writexl::write_xlsx(vals, path = filename_omega)
        }
      }
    }
    
    dev.off()
    
  } else {
    warning(paste("Model file not found for run:", run_name, "- Skipping."))
  }
}