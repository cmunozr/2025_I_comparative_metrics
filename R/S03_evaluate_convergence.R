# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Fitted models
#	OUTPUT. MCMC convergence statistics for selected model parameters,
# illustrated (for all RUNs performed thus far in S3) in the file "results/MCMC_convergence.pdf",
# and the text file "results/MCMC_convergence.txt".
# Also generates trace plots for selected parameters saved in "results/MCMC_trace_plots.pdf".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
showBeta = NULL #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = NULL #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = NULL #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = NULL #Default: convergence of Omega shown for 50 randomly selected species pairs
showRho = NULL #Default: showRho = FALSE, convergence not shown for rho-parameters
showAlpha = NULL #Default: showAlpha = FALSE, convergence not shown for alpha-parameters
showV = NULL #Default: showV = FALSE, convergence not shown for V-parameters
num_trace_plots = NULL #Default: num_trace_plots = 5, number of random trace plots per parameter group
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100 # Subsample 100 elements if Omega matrix is larger than 100*100
showRho = TRUE
showAlpha = TRUE
showV = TRUE
num_trace_plots = 5 # Number of random parameters to plot trace for in each group
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
modelFile = "fitted_fbsF_001.RData"
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(showBeta)) showBeta = TRUE
if(is.null(showGamma)) showGamma = FALSE
if(is.null(showOmega)) showOmega = FALSE
if(is.null(maxOmega)) maxOmega = 50
if(is.null(showRho)) showRho = FALSE
if(is.null(showAlpha)) showAlpha = FALSE
if(is.null(showV)) showV = FALSE
if(is.null(num_trace_plots)) num_trace_plots = 5


library(Hmsc)
library(colorspace)
library(vioplot)
library(coda) # Explicitly load coda for clarity and ESS

# Define samples and thin
samples_list = 1000
thin_list = 100
nst = length(thin_list)
nChains = 4

text.file = file.path(resultDir, "/MCMC_convergence.txt")
pdf.file.psrf = file.path(resultDir, "/MCMC_convergence.pdf")
pdf.file.ess = file.path(resultDir, "/MCMC_ESS.pdf")
pdf.file.trace = file.path(resultDir, "/MCMC_trace_plots.pdf")

cat("MCMC Convergence Statistics\n\n", file=text.file, sep="")

# Initialize lists to store PSRF and ESS for plotting
psrf.beta.all = list()
ess.beta.all = list()
psrf.gamma.all = list()
ess.gamma.all = list()
psrf.omega.all = list()
ess.omega.all = list()
psrf.alpha.all = list()
ess.alpha.all = list()
psrf.rho.all = list()
ess.rho.all = list()
psrf.V.all = list()
ess.V.all = list()

# Open PDF for trace plots
pdf(file = pdf.file.trace)

Lst = 1
while(Lst <= nst){ # Loops through thinning/sample settings
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(modelDir, modelFile)
  if(file.exists(filename)){
    load(filename)
    cat(c("\n", filename, "\n\n"), file=text.file, sep="", append=TRUE)
    nm = length(models)
    
    for(j in 1:nm){ 
      model_name = names(models)[j]
      cat(c("\n--- Model: ", model_name, " ---\n\n"), file=text.file, sep="", append=TRUE)
      
      # Convert to coda object
      # Use spNamesNumbers and covNamesNumbers for more detailed output in coda objects
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(TRUE,FALSE), covNamesNumbers = c(TRUE,FALSE))
      nr = models[[j]]$nr # Number of random levels
      
      cat(paste0("MCMC parameters: Thin=", thin, ", Samples=", samples, ", Chains=", nChains, "\n\n"), file=text.file, sep="", append=TRUE)
      
      sum_text = "Min.,   1st Qu.   Median     Mean     3rd Qu.    Max"
      
      # --- BETA ---
      if(showBeta){
        cat("\n--- Beta (Species Niches) ---\n", file=text.file, sep="", append=TRUE)
        psrf = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
        ess = effectiveSize(mpost$Beta)
        
        cat("\nPSRF Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(psrf[,1]), file=text.file, append=TRUE)
        cat("\nESS Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(ess), file=text.file, append=TRUE)
        
        psrf.beta.all[[model_name]] = psrf[,1]
        ess.beta.all[[model_name]] = ess
        
        # Selective Trace Plots for Beta
        par_names_beta = colnames(mpost$Beta[[1]])
        if(length(par_names_beta) > 0 && num_trace_plots > 0) {
          num_to_plot = min(num_trace_plots, length(par_names_beta))
          selected_indices = sample(1:length(par_names_beta), num_to_plot, replace = FALSE)
          for (i in 1:length(selected_indices)) {
            p_index = selected_indices[i]
            p_name = par_names_beta[p_index] # Get the full parameter name for reference if needed
            
            main_title <- paste0("Trace of Beta: Parameter ", p_index) # Shorter main title
            
            plot(mpost$Beta[, p_name, drop = FALSE], main = main_title)
            mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
          }
        }
      }
      
      # --- GAMMA ---
      if(showGamma){
        cat("\n--- Gamma (Influence of Traits on Niches) ---\n", file=text.file, sep="", append=TRUE)
        psrf = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
        ess = effectiveSize(mpost$Gamma)
        
        cat("\nPSRF Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(psrf[,1]), file=text.file, append=TRUE)
        cat("\nESS Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(ess), file=text.file, append=TRUE)
        
        psrf.gamma.all[[model_name]] = psrf[,1]
        ess.gamma.all[[model_name]] = ess
        
        # Selective Trace Plots for Gamma
        par_names_gamma = colnames(mpost$Gamma[[1]])
        if(length(par_names_gamma) > 0 && num_trace_plots > 0) {
          num_to_plot = min(num_trace_plots, length(par_names_gamma))
          selected_indices = sample(1:length(par_names_gamma), num_to_plot, replace = FALSE)
          for (i in 1:length(selected_indices)) {
            p_index = selected_indices[i]
            p_name = par_names_gamma[p_index]
            
            main_title <- paste0("Trace of Gamma: Parameter ", p_index)
            
            plot(mpost$Gamma[, p_name, drop = FALSE], main = main_title)
            mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
          }
        }
      }
      
      # --- RHO ---
      if(showRho & !is.null(mpost$Rho)){
        cat("\n--- Rho (Phylogenetic Signal) ---\n", file=text.file, sep="", append=TRUE)
        psrf = gelman.diag(mpost$Rho, multivariate=FALSE)$psrf
        ess = effectiveSize(mpost$Rho)
        
        cat("\nPSRF:\n", file=text.file, append=TRUE)
        cat(psrf[1], file=text.file, sep="\n", append=TRUE)
        cat("\nESS:\n", file=text.file, append=TRUE)
        cat(ess[1], file=text.file, sep="\n", append=TRUE)
        
        psrf.rho.all[[model_name]] = psrf[1]
        ess.rho.all[[model_name]] = ess[1]
        
        # Selective Trace Plots for Rho
        if(num_trace_plots > 0) {
          main_title <- "Trace of Rho"
          plot(mpost$Rho, main = main_title)
          mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
        }
      }
      
      # --- V ---
      if(showV & !is.null(mpost$V)){
        cat("\n--- V (Residual Covariance of Species Niches) ---\n", file=text.file, sep="", append=TRUE)
        psrf = gelman.diag(mpost$V, multivariate=FALSE)$psrf
        ess = effectiveSize(mpost$V)
        
        cat("\nPSRF Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(psrf[,1]), file=text.file, append=TRUE)
        cat("\nESS Summary:\n", file=text.file, append=TRUE)
        cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
        cat(summary(ess), file=text.file, append=TRUE)
        
        psrf.V.all[[model_name]] = psrf[,1]
        ess.V.all[[model_name]] = ess
        
        # Selective Trace Plots for V (flatting to vectors for plotting)
        par_names_V = colnames(mpost$V[[1]])
        if(length(par_names_V) > 0 && num_trace_plots > 0) {
          num_to_plot = min(num_trace_plots, length(par_names_V))
          selected_indices = sample(1:length(par_names_V), num_to_plot, replace = FALSE)
          selected_par_names = par_names_V[selected_indices]
          for (i in 1:length(selected_indices)) {
            p_index = selected_indices[i]
            p_name = par_names_V[p_index]
            
            main_title <- paste0("Trace of V: Parameter ", p_index)
            
            plot(mpost$V[, p_name, drop = FALSE], main = main_title)
            mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
          }
        }
      }
      
      
      # --- OMEGA ---
      if(showOmega & nr > 0){
        cat("\n--- Omega (Species Associations) ---\n", file=text.file, sep="", append=TRUE)
        for(k in 1:nr){ # Loops through each random level
          random_level_name = names(models[[j]]$ranLevels)[k]
          cat(c("\nRandom Level: ", random_level_name, "\n"), file=text.file, sep="", append=TRUE)
          
          tmp = mpost$Omega[[k]] # Omega for the current random level
          z = ncol(tmp[[1]]) # Number of pairs (ns*ns)
          
          tmp_subsampled = tmp # Default to full
          if(z > maxOmega * maxOmega){ # Check against maxOmega * maxOmega because maxOmega refers to elements
            selected_col_indices = sample(1:z, size = maxOmega * maxOmega, replace = FALSE)
            tmp_subsampled = lapply(tmp, function(chain_mat) {
              chain_mat[, selected_col_indices, drop = FALSE]
            })
            class(tmp_subsampled) <- "mcmc.list"
          }
          
          psrf = gelman.diag(tmp_subsampled, multivariate = FALSE)$psrf
          ess = effectiveSize(tmp_subsampled)
          
          cat("\nPSRF Summary:\n", file=text.file, append=TRUE)
          cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
          cat(summary(psrf[,1]), file=text.file, append=TRUE)
          cat("\nESS Summary:\n", file=text.file, append=TRUE)
          cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
          cat(summary(ess), file=text.file, append=TRUE)
          
          # Store for plotting, adding random level name to the list name
          psrf.omega.all[[paste0(model_name, "_", random_level_name)]] = psrf[,1]
          ess.omega.all[[paste0(model_name, "_", random_level_name)]] = ess
          
          # Selective Trace Plots for Omega (subsample individual pairs)
          par_names_omega = colnames(mpost$Omega[[k]][[1]]) # Names are like "sp1:sp1", "sp1:sp2"
          if(length(par_names_omega) > 0 && num_trace_plots > 0) {
            num_to_plot = min(num_trace_plots, length(par_names_omega))
            selected_indices = sample(1:length(par_names_omega), num_to_plot, replace = FALSE)
            for (i in 1:length(selected_indices)) {
              p_index = selected_indices[i]
              p_name = par_names_omega[p_index] # Full original name like "sp1:sp2"
              
              main_title <- paste0("Trace of Omega: Pair ", p_index, " (", random_level_name, ")")
              
              plot(mpost$Omega[[k]][, p_name, drop = FALSE], main = main_title)
              mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
            }
          }
        }
      }
      
      # --- ALPHA ---
      if(showAlpha & nr > 0){
        cat("\n--- Alpha (Spatial Scale) ---\n", file=text.file, sep="", append=TRUE)
        for(k in 1:nr){ # Loops through each random level
          if(models[[j]]$ranLevels[[k]]$sDim > 0){ # Check if it's a spatial random level
            random_level_name = names(models[[j]]$ranLevels)[k]
            cat(c("\nRandom Level: ", random_level_name, "\n"), file=text.file, sep="", append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]], multivariate = FALSE)$psrf
            ess = effectiveSize(mpost$Alpha[[k]])
            
            cat("\nPSRF Summary:\n", file=text.file, append=TRUE)
            cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
            cat(summary(psrf[,1]), file=text.file, append=TRUE)
            cat("\nESS Summary:\n", file=text.file, append=TRUE)
            cat(paste(sum_text, "\n"), file=text.file, append=TRUE)
            cat(summary(ess), file=text.file, append=TRUE)
            
            # Store for plotting, adding random level name
            psrf.alpha.all[[paste0(model_name, "_", random_level_name)]] = psrf[,1]
            ess.alpha.all[[paste0(model_name, "_", random_level_name)]] = ess
            
            # Selective Trace Plots for Alpha
            par_names_alpha = colnames(mpost$Alpha[[k]][[1]])
            if(length(par_names_alpha) > 0 && num_trace_plots > 0) {
              num_to_plot = min(num_trace_plots, length(par_names_alpha))
              selected_indices = sample(1:length(par_names_alpha), num_to_plot, replace = FALSE)
              for (i in 1:length(selected_indices)) {
                p_index = selected_indices[i]
                p_name = par_names_alpha[p_index] # Full original name like "factor1"
                
                main_title <- paste0("Trace of Alpha: Factor ", p_index, " (", random_level_name, ")")
                
                plot(mpost$Alpha[[k]][, p_name, drop = FALSE], main = main_title)
                mtext(model_name, side = 3, line = 1, cex = 0.8) # Subtitle: model name
              }
            }
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

dev.off() # Close trace plot PDF

# --- Generate combined PSRF and ESS plots ---
# This part uses the collected lists to make violin plots
model_colors = rainbow_hcl(length(names(models)))
names(model_colors) = names(models)

# Open PDF for PSRF plots
pdf(file = pdf.file.psrf)

plot_vioplot_section <- function(data_list, title_suffix, y_range_full, y_range_zoom, color_map) {
  if(length(data_list) > 0) {
    # Convert list of vectors to a matrix for vioplot
    max_len = max(sapply(data_list, length))
    padded_data = lapply(data_list, function(x) c(x, rep(NA, max_len - length(x))))
    plot_matrix = do.call(cbind, padded_data)
    colnames(plot_matrix) = names(data_list)
    
    # Get colors for each column
    column_colors = sapply(names(data_list), function(n) {
      model_base_name = strsplit(n, "_")[[1]][1] # Extract base name if random level suffix is added
      if (model_base_name %in% names(color_map)) {
        return(color_map[model_base_name])
      } else {
        for (mn in names(color_map)) {
          if (grepl(mn, n)) {
            return(color_map[mn])
          }
        }
      }
      return("grey") # Default if no match
    })
    
    par(mfrow=c(2,1))
    vioplot(plot_matrix, col=column_colors, names=names(data_list), ylim=y_range_full, main=paste0("PSRF (", title_suffix, ")"))
    legend("topright", legend = names(models), fill=model_colors)
    vioplot(plot_matrix, col=column_colors, names=names(data_list), ylim=y_range_zoom, main=paste0("PSRF (", title_suffix, ") (Zoomed)"))
    legend("topright", legend = names(models), fill=model_colors)
  }
}

plot_vioplot_section_ess <- function(data_list, title_suffix, y_range_full, color_map, min_ess_threshold = 1000) { # Adjusted threshold for ESS
  if(length(data_list) > 0) {
    max_len = max(sapply(data_list, length))
    padded_data = lapply(data_list, function(x) c(x, rep(NA, max_len - length(x))))
    plot_matrix = do.call(cbind, padded_data)
    colnames(plot_matrix) = names(data_list)
    
    column_colors = sapply(names(data_list), function(n) {
      model_base_name = strsplit(n, "_")[[1]][1] # Extract base name if random level suffix is added
      if (model_base_name %in% names(color_map)) {
        return(color_map[model_base_name])
      } else {
        for (mn in names(color_map)) {
          if (grepl(mn, n)) {
            return(color_map[mn])
          }
        }
      }
      return("grey")
    })
    
    par(mfrow=c(1,1))
    vioplot(plot_matrix, col=column_colors, names=names(data_list), ylim=y_range_full, main=paste0("ESS (", title_suffix, ")"))
    abline(h = min_ess_threshold, col = "red", lty = 2) # Add a line for the ESS threshold
    text(x = par("usr")[2], y = min_ess_threshold, labels = paste("Min ESS:", min_ess_threshold), pos = 2, col = "red")
    legend("topright", legend = names(models), fill=model_colors)
  }
}


# Plotting PSRF
plot_vioplot_section(psrf.beta.all, "Beta", c(0, max(unlist(psrf.beta.all), na.rm = TRUE)), c(0.9,1.1), model_colors)
plot_vioplot_section(psrf.gamma.all, "Gamma", c(0, max(unlist(psrf.gamma.all), na.rm = TRUE)), c(0.9,1.1), model_colors)
plot_vioplot_section(psrf.V.all, "V", c(0, max(unlist(psrf.V.all), na.rm = TRUE)), c(0.9,1.1), model_colors)
plot_vioplot_section(psrf.omega.all, "Omega", c(0, max(unlist(psrf.omega.all), na.rm = TRUE)), c(0.9,1.1), model_colors)
#plot_vioplot_section(psrf.alpha.all, "Alpha", c(0, max(unlist(psrf.alpha.all), na.rm = TRUE)), c(0.9,1.1), model_colors)
# Rho is typically a single value, so a violin plot might not be the best representation, but keeping for consistency
if(length(psrf.rho.all) > 0) {
  # For single values, convert to a single-column matrix for vioplot
  rho_psrf_data = do.call(cbind, psrf.rho.all)
  colnames(rho_psrf_data) = names(psrf.rho.all)
  par(mfrow=c(2,1))
  vioplot(rho_psrf_data, col=model_colors[names(psrf.rho.all)], names=names(psrf.rho.all), ylim=c(0,max(unlist(psrf.rho.all), na.rm=TRUE)), main="PSRF (Rho)")
  legend("topright", legend = names(models), fill=model_colors)
  vioplot(rho_psrf_data, col=model_colors[names(psrf.rho.all)], names=names(psrf.rho.all), ylim=c(0.9,1.1), main="PSRF (Rho) (Zoomed)")
  legend("topright", legend = names(models), fill=model_colors)
}
dev.off() 

# Open PDF for ESS plots
pdf(file = pdf.file.ess)
plot_vioplot_section_ess(ess.beta.all, "Beta", c(0, max(unlist(ess.beta.all), na.rm = TRUE)), model_colors)
plot_vioplot_section_ess(ess.gamma.all, "Gamma", c(0, max(unlist(ess.gamma.all), na.rm = TRUE)), model_colors)
plot_vioplot_section_ess(ess.V.all, "V", c(0, max(unlist(ess.V.all), na.rm = TRUE)), model_colors)
plot_vioplot_section_ess(ess.omega.all, "Omega", c(0, max(unlist(ess.omega.all), na.rm = TRUE)), model_colors)
plot_vioplot_section_ess(ess.alpha.all, "Alpha", c(0, max(unlist(ess.alpha.all), na.rm = TRUE)), model_colors)
if(length(ess.rho.all) > 0) {
  rho_ess_data = do.call(cbind, ess.rho.all)
  colnames(rho_ess_data) = names(ess.rho.all)
  par(mfrow=c(1,1))
  vioplot(rho_ess_data, col=model_colors[names(ess.rho.all)], names=names(ess.rho.all), ylim=c(0,max(unlist(ess.rho.all), na.rm=TRUE)), main="ESS (Rho)")
  abline(h = 1000, col = "red", lty = 2) # Add a line for the ESS threshold
  text(x = par("usr")[2], y = 1000, labels = "Min ESS: 1000", pos = 2, col = "red")
  legend("topright", legend = names(models), fill=model_colors)
}
dev.off()
