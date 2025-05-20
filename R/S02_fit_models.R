nChains <- 4
nParallel <- 4 #Default: nParallel = nChains, 1 deactivated
localDir <-  "C:/Users/Carlos Munoz/Documents/Ph.D/6_courses/2025_I_comparative_metrics"
modelDir <-  file.path(localDir, "models")
dir.create(modelDir, recursive = TRUE, showWarnings = F)

library(Hmsc)
library(doParallel)

unfitted_models_file <- file.path(modelDir, "unfitted_models.RData")
load(file = unfitted_models_file)

nm <- 1#length(models)

# --- Define MCMC Sampling Parameters ---

if (.Platform$OS.type == "windows") {
  samples_list <-  c(1000) #, 250, 250) 
  thin_list <- c(1000) #, 1, 10)
} else {
  samples_list <- c(5,1000)  # 100,1000, 10000)
  thin_list <- c(1, 1000) # 100,1000, 10000)
}

# --- Setup Parallel Processing ---¨

#if (is.null(nParallel)) nParallel <- nChains
#registerDoParallel(cores = nParallel)

# --- Main Loop for Fitting Models ---

Lst = 1
while (Lst <= length(samples_list)) {
  thin <-  thin_list[Lst]
  samples <- samples_list[Lst]
  
  print(paste0("--- Starting: thin = ", thin, ", samples = ", samples, " ---"))
  
  for (mi in 1:nm) {
    # mi <- 1
    model_name <- names(models)[mi]
    # Define filename for this specific model and parameters
    filename <-  file.path(modelDir, paste0(
      model_name,
      "_thin_", as.character(thin),
      "_samples_", as.character(samples),
      "_chains_", as.character(nChains),
      ".Rdata"
    ))
    
    if (file.exists(filename)) {
      print(paste("Skipping (exists):", basename(filename)))
    } else {
      
      print(paste("Processing model:", model_name, " | Saving to:", basename(filename)))
      start_time <- Sys.time()
      print(start_time)
      
      m_unfitted <- models[[mi]]
      
      # Fit the model using MCMC with error handling
      m_fitted <- tryCatch({
        sampleMcmc(m_unfitted,
                   samples = samples,
                   thin = thin,
                   adaptNf = rep(ceiling(0.4 * samples * thin), m_unfitted$nr),
                   transient = ceiling(0.5 * samples * thin),
                   nChains = nChains,
                   nParallel = nParallel,
                   verbose = ceiling(samples / 10) # Keep some progress output
        )
      }, error = function(e) {
        print(paste("ERROR sampling model:", model_name))
        print(e$message) # Print just the error message
        NULL # Return NULL on error
      })
      
      # Save the individual fitted model if sampling was successful
      if (!is.null(m_fitted)) {
        save(m_fitted, file = filename)
        end_time <- Sys.time()
        print(paste("Finished model:", model_name, "| Time:", format(end_time - start_time)))
      } else {
        print(paste("Failed model:", model_name, "- File not saved."))
      }
    }
  }
  
  Lst <- Lst + 1
}

# --- Cleanup ---
print("--- Script finished ---")
