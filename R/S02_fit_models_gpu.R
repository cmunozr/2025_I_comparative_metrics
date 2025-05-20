# install.packages("devtools") # if not yet installed
# library(devtools)
# install_github("hmsc-r/HMSC")

library(Hmsc)
library(jsonify)

nChains <- 4
nParallel <- 4 #Default: nParallel = nChains, 1 deactivated
localDir <-  "C:/Users/Carlos Munoz/Documents/Ph.D/6_courses/2025_I_comparative_metrics"
modelDir <-  file.path(localDir, "models")
dir.create(modelDir, recursive = TRUE, showWarnings = F)

output_cmds_file <- file.path(modelDir, "python_commands_log.txt")


unfitted_models_file <- file.path(modelDir, "unfitted_models_small_NGPP.RData")
load(file = unfitted_models_file)

nm <- 1#length(models)

# --- Define MCMC Sampling Parameters ---

if (.Platform$OS.type == "windows") {
  samples_list <-  c(5, 250, 250, 1000, 1000) #, 250, 250) 
  thin_list <- c(1, 1, 250, 100, 1000) #, 1, 10)
} else {
  samples_list <- c(5)  # 100,1000, 10000)
  thin_list <- c(1) # 100,1000, 10000)
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
    model_sample <- paste0(
      model_name, "_GPU",
      "_thin_", as.character(thin),
      "_samples_", as.character(samples),
      "_chains_", as.character(nChains)
    )
    filename <-  file.path(modelDir, model_sample)
    
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
                   verbose = ceiling(samples / 10), # Keep some progress output
                   engine="HPC",
                   updater=list(GammaEta=FALSE)
        )
      }, error = function(e) {
        print(paste("ERROR sampling model:", model_name))
        print(e$message) # Print just the error message
        NULL # Return NULL on error
      })
      
      # Save the individual fitted model if sampling was successful
      if (!is.null(m_fitted)) {
        saveRDS(to_json(m_fitted), file=paste0(filename, ".rds"))
        end_time <- Sys.time()
        print(paste("Finished model:", model_name, "| Time:", format(end_time - start_time)))
        
        post_file_path <-  file.path(paste0("models/", model_sample, "_post_file.rds"))
        python_cmd_args <- paste("-m hmsc.run_gibbs_sampler",
                                 "--input", shQuote(paste0("models/", model_sample, ".rds")),
                                 "--output", shQuote(post_file_path),
                                 "--samples", samples,
                                 "--transient", ceiling(0.5 * samples * thin),
                                 "--thin", thin,
                                 "--verbose", ceiling(samples / 10),
                                 "--fp 32", 
                                 "&> nohupasimov.out &")
        full_python_command <- paste("nohup python", python_cmd_args)
        
        # Guardar el comando en el archivo de texto
        # 'append = TRUE' para añadir a las líneas existentes
        # 'fill = TRUE' o añadir '\n' para asegurar nueva línea por comando
        cat(full_python_command, file = output_cmds_file, append = TRUE, sep = "\n")
        
        
        
      } else {
        print(paste("Failed model:", model_name, "- File not saved."))
      }
    }
  }
  
  Lst <- Lst + 1
}



# --- Cleanup ---
print("--- Script finished ---")
