# install.packages("devtools") # if not yet installed
# library(devtools)
# install_github("hmsc-r/HMSC")

library(Hmsc)
library(jsonify)

nChains <- 4
nParallel <- 4 #Default: nParallel = nChains, 1 deactivated
gpu_Parallel <- TRUE ### <<<<<---------------------------------
localDir <- "C:/Users/Carlos Munoz/Documents/Ph.D/6_courses/2025_I_comparative_metrics"
modelDir <- file.path(localDir, "models")
dir.create(modelDir, recursive = TRUE, showWarnings = F)

name_unfitted_models <- "unfitted_fbsF_001"
unfitted_models_file <- file.path(modelDir, paste0(name_unfitted_models, ".RData"))
output_cmds_file <- file.path(modelDir, paste0("python_commands_",name_unfitted_models, "log.txt"))

load(file = unfitted_models_file)

nm <- length(models)

# --- Define MCMC Sampling Parameters ---

if (.Platform$OS.type == "windows") {
  samples_list <- c(5, 1000, 1000) #, 250, 250)
  thin_list <- c(1, 100, 1000) #, 1, 10)
} else {
  samples_list <- c(5) # 100,1000, 10000)
  thin_list <- c(1) # 100,1000, 10000)
}

# --- Setup Parallel Processing ---¨

#if (is.null(nParallel)) nParallel <- nChains
#registerDoParallel(cores = nParallel)

# --- Main Loop for Fitting Models ---

Lst = 1
while (Lst <= length(samples_list)) {
  thin <- thin_list[Lst]
  samples <- samples_list[Lst]
  
  print(paste0("--- Starting: thin = ", thin, ", samples = ", samples, " ---"))
  
  for (mi in 1:nm) {
    # mi <- 1
    model_name <- names(models)[mi]
    
    # Define the base name for the model's outputs, incorporating sampling parameters
    model_output_base <- paste0(
      model_name, "_GPU",
      "_thin_", as.character(thin),
      "_samples_", as.character(samples),
      "_chains_", as.character(nChains)
    )
    
    # Define the subdirectory path for this specific model and parameters
    # This creates a structure like: models/MyModel_GPU_thin_100_samples_1000_chains_4/
    current_model_subdir <- file.path(modelDir, model_output_base)
    dir.create(current_model_subdir, recursive = TRUE, showWarnings = F) # Create the subdirectory
    
    # Define filename for the unfitted model within its new subdirectory
    filename_unfitted_rds <- file.path(current_model_subdir, paste0(model_output_base, ".rds"))
    
    if (file.exists(filename_unfitted_rds)) {
      print(paste("Skipping (exists):", basename(filename_unfitted_rds)))
    } else {
      
      print(paste("Processing model:", model_name, " | Saving to:", current_model_subdir))
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
                   engine = "HPC",
                   updater = list(GammaEta = FALSE)
        )
      }, error = function(e) {
        print(paste("ERROR sampling model:", model_name))
        print(e$message) # Print just the error message
        NULL # Return NULL on error
      })
      
      # Save the individual fitted model if sampling was successful
      if (!is.null(m_fitted)) {
        
        saveRDS(to_json(m_fitted), file = filename_unfitted_rds)
        end_time <- Sys.time()
        print(paste("Finished model:", model_name, "| Time:", format(end_time - start_time)))
        
        if(!gpu_Parallel){
          post_file_path <-  file.path(paste0("models/", model_output_base, "_post_file.rds"))
          python_cmd_args <- paste("-m hmsc.run_gibbs_sampler",
                                   "--input", shQuote(paste0("models/", model_output_base, ".rds")),
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
          
        }else{
          
            chain_cmd_args_list <- vector("list", nChains)
            for (cInd in 1:nChains) {
              # Define a unique output file path for each chain's posterior within the subdirectory
              # This creates names like:
              # models/MyModel_GPU_thin_X_samples_Y_chains_Z/MyModel_GPU_thin_X_samples_Y_chains_Z_post_chain00_file.rds
              chain_post_file_path <- file.path(current_model_subdir, paste0(model_output_base, "_post_chain", sprintf("%.2d", cInd - 1), "_file.rds"))
              
              # Construct the Python command arguments
              input_path_for_python <- file.path("/scratch", "asimov", "munozcs", "models", model_output_base, paste0(model_output_base, ".rds")) #### MANUAL CHANGE
              output_path_for_python <- file.path("/scratch", "asimov", "munozcs", "models", model_output_base, paste0(model_output_base, "_post_chain", sprintf("%.2d", cInd - 1), "_file.rds")) #### MANUAL CHANGE
              
              input_path_for_python <- file.path("models", model_output_base, paste0(model_output_base, ".rds")) #### MANUAL CHANGE
              output_path_for_python <- file.path("models", model_output_base, paste0(model_output_base, "_post_chain", sprintf("%.2d", cInd - 1), "_file.rds")) #### MANUAL CHANGE
              
              chain_cmd_args <- paste("-m hmsc.run_gibbs_sampler",
                                      "--input", shQuote(input_path_for_python), # Initialized model from R
                                      "--output", shQuote(output_path_for_python),
                                      "--samples", samples,
                                      "--transient", ceiling(0.5 * samples * thin),
                                      "--thin", thin,
                                      "--verbose", ceiling(samples / 10),
                                      "--chain", cInd - 1, # Python uses 0-indexed chains
                                      "--fp 64"
              )
              
              # For the nohup log file path in the command itself, it also needs to be relative
              # to where the nohup command is executed or an absolute path.
              nohup_log_file_for_cmd <- file.path(paste0(model_output_base, "_chain", sprintf("%.2d", cInd - 1), "_nohup.out")) #### MANUAL CHANGE
              
              
              # Combine nohup, python executable, arguments, and output redirection
              full_python_command_with_nohup <- paste(
                "nohup python",
                chain_cmd_args,
                "&>",
                shQuote(nohup_log_file_for_cmd),
                "&"
              )
              
              # Store the command (e.g., to write to a log file or execute later)
              chain_cmd_args_list[[cInd]] <- full_python_command_with_nohup
              
              # --- Execute the command for each chain ---
              # Note: If we want to execute these commands directly from R, we would use system()
              # For example: system(full_python_command_with_nohup, wait = FALSE)
              # However, current setup writes them to a log file, which is good for HPC.
              cat(full_python_command_with_nohup, file = output_cmds_file, append = TRUE, sep = "\n")
              
          }
        
        }
      } else {
        print(paste("Failed model:", model_name, "- File not saved."))
      }
    }
  }
  
  Lst <- Lst + 1
}

# --- Cleanup ---
print("--- Script finished ---")