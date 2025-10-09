# --- 1. Load Required Libraries ---
library(Hmsc)
library(jsonify)

source("code/config_model.R")

# --- 2. Configuration and Setup ---

run_name <- generate_run_name(run_config)

# -- Execution Mode --
# 1: Simple TXT file (one command runs all chains).
# 2: Detailed TXT file (one command per chain).
# 3: Generate executable .sh scripts for each GPU.
execution_mode <- run_config$gpu$execution_mode 

# Define general paths
paths <- list(
  local_dir = getwd(),
  models_dir = file.path(getwd(), "models"),
  unfitted_models_file = file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
)

dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 4. Load the Unfitted Model ---
load(file = paths$unfitted_models_file)
n_models <- length(models) # Should typically be 1 if you define one model per file

# --- 5. Define the Parameter Grid for MCMC Sampling ---
mcmc_parameter_grid <- run_config$mcmc

# -----------------------------------------------
# Helper Function

generate_python_command <- function(input_path, output_path, params, chain_index = NULL, include_chain_arg = TRUE) {
  
  # Calculate numeric values first
  transient_value <- ceiling(params$transient_proportion * params$samples * params$thin)
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

# Initialize lists to store the generated commands
gpu_command_list <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)
all_commands_list <- list() # For modes 1 and 2

for (i in 1:nrow(mcmc_parameter_grid)) {
  # i <- 1
  run_name_i <- run_name[i]
  
  # define specific paths for each run
  paths_i <- list(
    run_specific_dir = file.path("models", run_name_i),
    unfitted_models_rds = file.path(file.path("models", run_name_i), paste0(run_name_i, ".rds")),
    python_commands_file = file.path(file.path("models", run_name_i), paste0("python_commands_", run_name_i, "log.txt"))
  )
  
  paths_i <- c(paths, paths_i)
  
  current_params <- mcmc_parameter_grid[i, ]
  message(paste0("\n--- Processing Parameter Set: thin = ", current_params$thin, ", samples = ", current_params$samples, " ---"))
  
  for (model_index in 1:n_models) {
    # model_index <- 1
    model_name <- run_config$model_id
    unfitted_model <- models[[model_index]]
    
    # --- A. Prepare and Save the .rds Model File ---
    model_output_base <- run_name_i
    current_model_subdir <- paths_i$run_specific_dir
    dir.create(current_model_subdir, recursive = TRUE, showWarnings = FALSE)
    filename_unfitted_rds <- paths_i$unfitted_models_rds
    
    model_prepared_successfully <- FALSE
    if (file.exists(filename_unfitted_rds)) {
      message(paste("Skipping .rds creation (already exists):", basename(filename_unfitted_rds)))
      model_prepared_successfully <- TRUE
    } else {
      message(paste("Processing and saving .rds for:", model_name))
      fitted_model <- tryCatch({
        sampleMcmc(
          unfitted_model, samples = current_params$samples, thin = current_params$thin,
          adaptNf = ceiling(current_params$adapt_nf_proportion * current_params$samples * current_params$thin),
          transient = ceiling(current_params$transient_proportion * current_params$samples * current_params$thin),
          nChains = current_params$n_chains, 
          nParallel = current_params$n_chains,
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
        input_path <- file.path(run_config$server$server_models_dir, run_name_i, basename(paths_i$unfitted_models_rds))
        output_path <- file.path(run_config$server$server_models_dir, run_name_i, paste0(model_output_base, "_post.rds"))
        
        python_command <- generate_python_command(input_path, output_path, params = current_params, include_chain_arg = FALSE)
        nohup_log_file <- paste0(model_output_base, "_nohup.out")
        full_shell_command <- paste("nohup", python_command, "&>", shQuote(nohup_log_file), "&")
        
        all_commands_list[[model_output_base]] <- full_shell_command
        
        # --- LOGIC FOR MODES 2 & 3 (This structure is now clean) ---
      } else if (execution_mode == 2 || execution_mode == 3) {
        if (execution_mode == 2) message("Mode 2: Generating detailed nohup commands...")
        if (execution_mode == 3) message("Mode 3: Generating commands for .sh scripts...")
        
        for (chain_index in 0:(current_params$n_chains - 1)) {
          # chain_index <- 1
          # Generate the core python command
          input_path <- file.path(run_config$server$server_models_dir, run_name_i, basename(paths_i$unfitted_models_rds))
          output_filename <- paste0(model_output_base, "_post_chain", sprintf("%.2d", chain_index), "_file.rds")
          output_path <- file.path(run_config$server$server_models_dir, run_name_i, output_filename)
          python_command <- generate_python_command(input_path, output_path, params = current_params, chain_index, include_chain_arg = TRUE)
          
          # Logic for Mode 2
          if (execution_mode == 2) {
            nohup_log_file <- paste0(model_output_base, "_chain", sprintf("%.2d", chain_index), "_nohup.out")
            full_shell_command <- paste("nohup", python_command, "&>", shQuote(nohup_log_file), "&")
            all_commands_list[[length(all_commands_list) + 1]] <- full_shell_command
          }
          
          # Logic for Mode 3
          if (execution_mode == 3) {
            nohup_log_name <- paste0(model_output_base, "_chain", sprintf("%.2d", chain_index), "_nohup.out")
            gpu_to_assign <- chain_index %% run_config$gpu$n_gpus_available
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
  con <- file(file.path(paths$models_dir, paste0(model_name, "_python_commands.txt")), open = "wb")
  cat(commands, file = con, sep = "\n")
  close(con)
  
  message(paste("Successfully created:", model_name, "python commands file"))
  
} else if (execution_mode == 3) {
    message("\n--- Generating .sh script for each GPU for Mode 3 ---")
    for (gpu_index in 0:(run_config$gpu$n_gpus_available - 1)) {
      sh_filename <- file.path(paths$models_dir, paste0(model_name, "_run_gpu_", gpu_index, ".sh"))
      commands_for_this_gpu <- gpu_command_list[[gpu_index + 1]]
      
      if (nrow(commands_for_this_gpu) == 0) { next }
      
      header <- c(
        "#!/bin/bash", "# Auto-generated by R script", "",
        paste0("PYTHON_ENV=\"", run_config$server$python_env, "\""),
        paste0("NOHUP_DIR=\"", run_config$server$nohup_dir, "\""),
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
