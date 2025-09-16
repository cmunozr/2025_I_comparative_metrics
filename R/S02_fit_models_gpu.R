# --- 1. Load Required Libraries ---
library(Hmsc)
library(jsonify)

# --- 2. Configuration and Setup ---

# -- Main Execution Mode --
# 1: Simple TXT file (one command runs all chains).
# 2: Detailed TXT file (one command per chain).
# 3: Generate executable .sh scripts for each GPU.
execution_mode <- 1

# -- HPC & MCMC Parameters --
config <- list(
  n_chains = 4,
  n_parallel = 4,
  n_gpus_available = 2,
  transient_proportion = 0.5,
  adapt_nf_proportion = 0.4
)

# -- Server Environment (for Mode 3) --
server_config <- list(
  python_env = "/home/avesta/munozcs/Documents/hmsc-venv",
  nohup_dir = "/home/avesta/munozcs/Documents/nohup_chains"
)

# -- File Paths --
local_dir <- getwd()
paths <- list(
  local_dir = local_dir,
  models_dir = file.path(local_dir, "models"),
  input_folder_server = "models",
  output_folder_server = "models"
)
dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# -- Model & Output File Definitions --
unfitted_models_filename <- "unfitted_fbsF_001"
paths$unfitted_models_file <- file.path(paths$models_dir, paste0(unfitted_models_filename, ".RData"))
paths$python_commands_file <- file.path(paths$models_dir, paste0("python_commands_", unfitted_models_filename, "log.txt"))

load(file = paths$unfitted_models_file)
n_models <- length(models)

# --- 3. Define the Parameter Grid for MCMC Sampling ---
mcmc_parameter_grid <- data.frame(
  samples = c(1000),
  thin = c(1000)
)

# -----------------------------------------------
# Helper Function

generate_python_command <- function(input_path, output_path, params, chain_index = NULL, include_chain_arg = TRUE) {
  
  # Calculate numeric values first
  transient_value <- ceiling(config$transient_proportion * params$samples * params$thin)
  verbose_value <- ceiling(params$samples / 10)
  
  python_args <- paste(
    "-m hmsc.run_gibbs_sampler", 
    "--input", shQuote(input_path), 
    "--output", shQuote(output_path),
    "--samples", format(params$samples, scientific = FALSE), 
    "--transient", format(transient_value, scientific = FALSE),
    "--thin", format(params$thin, scientific = FALSE), 
    "--verbose", format(verbose_value, scientific = FALSE), 
    
    "--fp 64"
  )
  
  if (include_chain_arg && !is.null(chain_index)) {
    python_args <- paste(python_args, "--chain", chain_index)
  }
  return(paste("python", python_args))
}
# -----------------------------------------------
# Main loop

# Initialize lists to store the generated commands
gpu_command_list <- rep(list(data.frame(command=character(), log_filename=character())), config$n_gpus_available)
all_commands_list <- list() # For modes 1 and 2

for (i in 1:nrow(mcmc_parameter_grid)) {
  current_params <- mcmc_parameter_grid[i, ]
  message(paste0("\n--- Processing Parameter Set: thin = ", current_params$thin, ", samples = ", current_params$samples, " ---"))
  
  for (model_index in 1:n_models) {
    model_name <- names(models)[model_index]
    unfitted_model <- models[[model_index]]
    
    # --- A. Prepare and Save the .rds Model File ---
    model_output_base <- paste0(model_name, "_GPU", "_thin_", current_params$thin, "_samples_", current_params$samples, "_chains_", config$n_chains)
    current_model_subdir <- file.path(paths$models_dir, model_output_base)
    dir.create(current_model_subdir, recursive = TRUE, showWarnings = FALSE)
    filename_unfitted_rds <- file.path(current_model_subdir, paste0(model_output_base, ".rds"))
    
    model_prepared_successfully <- FALSE
    if (file.exists(filename_unfitted_rds)) {
      message(paste("Skipping .rds creation (already exists):", basename(filename_unfitted_rds)))
      model_prepared_successfully <- TRUE
    } else {
      message(paste("Processing and saving .rds for:", model_name))
      fitted_model <- tryCatch({
        sampleMcmc(
          unfitted_model, samples = current_params$samples, thin = current_params$thin,
          adaptNf = rep(ceiling(config$adapt_nf_proportion * current_params$samples * current_params$thin), unfitted_model$nr),
          transient = ceiling(config$transient_proportion * current_params$samples * current_params$thin),
          nChains = config$n_chains, nParallel = config$n_parallel,
          verbose = ceiling(current_params$samples / 10), engine = "HPC", updater = list(GammaEta = FALSE)
        )
      }, error = function(e) {
        warning(paste("ERROR preparing model:", model_name, "\nMessage:", e$message))
        return(NULL)
      })
      
      if (!is.null(fitted_model)) {
        saveRDS(to_json(fitted_model), file = filename_unfitted_rds)
        message("Successfully saved .rds file.")
        model_prepared_successfully <- TRUE
      }
    }
    
    # --- B. Generate Commands Based on the Chosen Mode ---
    if (model_prepared_successfully) {
      
      # --- LOGIC FOR MODE 1 ---
      if (execution_mode == 1) {
        message("Mode 1: Generating simple nohup command...")
        input_path <- file.path(paths$input_folder_server, model_output_base, paste0(model_output_base, ".rds"))
        output_path <- file.path(paths$output_folder_server, model_output_base, paste0(model_output_base, "_post.rds"))
        
        python_command <- generate_python_command(input_path, output_path, current_params, include_chain_arg = FALSE)
        nohup_log_file <- paste0(model_output_base, "_nohup.out")
        full_shell_command <- paste("nohup", python_command, "&>", shQuote(nohup_log_file), "&")
        
        all_commands_list[[model_output_base]] <- full_shell_command
        
        # --- LOGIC FOR MODES 2 & 3 (This structure is now clean) ---
      } else if (execution_mode == 2 || execution_mode == 3) {
        if (execution_mode == 2) message("Mode 2: Generating detailed nohup commands...")
        if (execution_mode == 3) message("Mode 3: Generating commands for .sh scripts...")
        
        # This is the single, correct loop for the chains. The extra nested loop was removed.
        for (chain_index in 0:(config$n_chains - 1)) {
          
          # Generate the core python command
          input_path <- file.path(paths$input_folder_server, model_output_base, paste0(model_output_base, ".rds"))
          output_filename <- paste0(model_output_base, "_post_chain", sprintf("%.2d", chain_index), "_file.rds")
          output_path <- file.path(paths$output_folder_server, model_output_base, output_filename)
          python_command <- generate_python_command(input_path, output_path, current_params, chain_index, include_chain_arg = TRUE)
          
          # Logic for Mode 2
          if (execution_mode == 2) {
            nohup_log_file <- paste0(model_output_base, "_chain", sprintf("%.2d", chain_index), "_nohup.out")
            full_shell_command <- paste("nohup", python_command, "&>", shQuote(nohup_log_file), "&")
            all_commands_list[[length(all_commands_list) + 1]] <- full_shell_command
          }
          
          # Logic for Mode 3
          if (execution_mode == 3) {
            nohup_log_name <- paste0(model_output_base, "_chain", sprintf("%.2d", chain_index), "_nohup.out")
            gpu_to_assign <- chain_index %% config$n_gpus_available
            new_entry <- data.frame(command = python_command, log_filename = nohup_log_name)
            gpu_command_list[[gpu_to_assign + 1]] <- rbind(gpu_command_list[[gpu_to_assign + 1]], new_entry)
          }
        }
      }
    }
  }
}

# -----------------------------------------------
# File generation

if (execution_mode == 1 || execution_mode == 2) {
  # For modes 1 & 2, write all commands to a single log file
  message(paste0("\n--- Generating master command log file for Mode ", execution_mode, " ---"))
  commands <- c(
    unlist(all_commands_list),
    "",
    "For manual run use:",
    "CUDA_VISIBLE_DEVICES=0"
  )
  con <- file(paths$python_commands_file, open = "wb")
  cat(commands, file = con, sep = "\n")
  close(con)
  
  message(paste("Successfully created:", paths$python_commands_file))
  
} else if (execution_mode == 3) {
    message("\n--- Generating .sh script for each GPU for Mode 3 ---")
    for (gpu_index in 0:(config$n_gpus_available - 1)) {
      sh_filename <- file.path(paths$models_dir, paste0("run_gpu_", gpu_index, ".sh"))
      commands_for_this_gpu <- gpu_command_list[[gpu_index + 1]]
      
      if (nrow(commands_for_this_gpu) == 0) { next }
      
      header <- c(
        "#!/bin/bash", "# Auto-generated by R script", "",
        paste0("PYTHON_ENV=\"", server_config$python_env, "\""),
        paste0("NOHUP_DIR=\"", server_config$nohup_dir, "\""),
        "mkdir -p \"$NOHUP_DIR\"", "",
        "echo \"Activating Python environment\"", paste0("source \"$PYTHON_ENV/bin/activate\""), "",
        paste0("export CUDA_VISIBLE_DEVICES=", gpu_index), "echo \"Running commands on GPU: $CUDA_VISIBLE_DEVICES\"", ""
      )
      
      wrapped_commands <- character(nrow(commands_for_this_gpu))
      for (j in 1:nrow(commands_for_this_gpu)) {
        cmd <- commands_for_this_gpu$command[j]
        log_name <- commands_for_this_gpu$log_filename[j]
        
        wrapped_commands[j] <- paste0(
          "echo \"----------------------------------------------------\"\n",
          "echo \"Starting: ", log_name, "\"\n",
          "nohup ", cmd, " &> \"$NOHUP_DIR/", log_name, "\" &\n",
          "pid=$!\n",
          "echo \"Process started with PID: $pid\"\n",
          "wait $pid\n",
          "echo \"Process $pid finished.\"\n"
        )
      }
      
      footer <- c("", paste0("echo \"All tasks for GPU ", gpu_index, " are complete.\""))
      
      # Write the file in binary mode to ensure correct line endings
      con <- file(sh_filename, open = "wb")
      cat(c(header, wrapped_commands, footer), file = con, sep = "\n")
      close(con)
    
    message(paste("Successfully created sh script:", sh_filename))
  }
}
