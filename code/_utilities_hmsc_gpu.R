# A set of functions to prepare and generate run scripts for Hmsc-HPC (GPU)

#' 1. Prepare the Hmsc Model Object for HPC
#'
#' This function takes an unfitted Hmsc model and the MCMC parameters,
#' and runs sampleMcmc with engine="HPC" to prepare the object for
#' parallel GPU execution.
#'
#' @param unfitted_model The loaded, unfitted Hmsc model object.
#' @param mcmc_params A list or data.frame row with MCMC settings
#'        (e.g., samples, thin, n_chains, transient_proportion, etc.)
#' @return A prepared Hmsc model object (or NULL on error).
prepare_hpc_model <- function(unfitted_model, mcmc_params) {
  
  # Calculate transient and verbose from proportions
  transient_val <- ceiling(mcmc_params$transient_proportion * mcmc_params$samples * mcmc_params$thin)
  adapt_nf_val <- ceiling(mcmc_params$adapt_nf_proportion * mcmc_params$samples * mcmc_params$thin)
  verbose_val <- ceiling(mcmc_params$samples / 10)
  
  prepared_model <- tryCatch({
    Hmsc::sampleMcmc(
      unfitted_model, 
      samples = mcmc_params$samples, 
      thin = mcmc_params$thin,
      adaptNf = adapt_nf_val,
      transient = transient_val,
      nChains = mcmc_params$n_chains, 
      nParallel = mcmc_params$n_chains,
      verbose = verbose_val, 
      engine = "HPC", 
      updater = list(GammaEta = FALSE)
    )
  }, error = function(e) {
    warning(paste("ERROR preparing model:", e$message))
    return(NULL)
  })
  
  return(prepared_model)
}


#' 2. Save the Prepared Model as a JSON-formatted .rds
#'
#' This function takes the prepared model from prepare_hpc_model()
#' and saves it to disk in the format required by the Python sampler.
#'
#' @param prepared_model The model object from prepare_hpc_model().
#' @param output_path The full file path to save the .rds file.
#' @return TRUE on success, FALSE on failure.
save_prepared_model <- function(prepared_model, output_path) {
  if (is.null(prepared_model)) {
    warning("Model is NULL. Cannot save.")
    return(FALSE)
  }
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(output_path)) {
    message("Skipping .rds creation (already exists): ", basename(output_path))
    return(TRUE)
  }
  
  message("Saving prepared model to: ", output_path)
  tryCatch({
    saveRDS(jsonify::to_json(prepared_model), file = output_path)
    message("Successfully saved .rds file.")
    return(TRUE)
  }, error = function(e) {
    warning(paste("ERROR saving .rds file:", e$message))
    return(FALSE)
  })
}


#' 3. Generate All Run Scripts/Commands
#'
#' Generates the command lists (for all modes) for a *single* model setup.
#' Does NOT write any files.
#'
#' @param base_model_name The unique name for THIS run (e.g., "m_sp_tr_co_all_S1000_T100_C4").
#' @param run_specific_dir_server The full path *on the server* where the .rds file is.
#' @param mcmc_params A list with MCMC settings (samples, thin, n_chains, etc.).
#' @param run_config The full config list.
#' @return A list containing $all_commands (for modes 1/2) and $gpu_commands (for mode 3).
generate_commands <- function(base_model_name, run_specific_dir_server, mcmc_params, run_config) {
  
  execution_mode <- run_config$gpu$execution_mode
  n_chains <- mcmc_params$n_chains
  
  # Initialize empty lists
  all_commands_list <- list() 
  gpu_command_list <- rep(list(data.frame(command=character(), log_filename=character())), run_config$gpu$n_gpus_available)
  
  model_output_base <- base_model_name
  input_rds_name <- paste0(base_model_name, ".rds")
  
  # --- LOGIC FOR MODE 1 ---
  if (execution_mode == 1) {
    message("Mode 1: Generating command for: ", base_model_name)
    input_path <- file.path(run_specific_dir_server, input_rds_name)
    output_path <- file.path(run_specific_dir_server, paste0(model_output_base, "_post.rds"))
    
    python_command <- .generate_python_command(input_path, output_path, params = mcmc_params, include_chain_arg = FALSE)
    nohup_log_file <- paste0(model_output_base, "_nohup.out")
    full_shell_command <- paste("nohup", python_command, "&>", shQuote(nohup_log_file), "&")
    
    all_commands_list[[model_output_base]] <- full_shell_command
    
    # --- LOGIC FOR MODES 2 & 3 ---
  } else if (execution_mode == 2 || execution_mode == 3) {
    if (execution_mode == 2) message("Mode 2: Generating commands for: ", base_model_name)
    if (execution_mode == 3) message("Mode 3: Generating commands for: ", base_model_name)
    
    for (chain_index in 0:(n_chains - 1)) {
      input_path <- file.path(run_specific_dir_server, input_rds_name)
      output_filename <- paste0(model_output_base, "_post_chain", sprintf("%.2d", chain_index), "_file.rds")
      output_path <- file.path(run_specific_dir_server, output_filename)
      python_command <- .generate_python_command(input_path, output_path, params = mcmc_params, chain_index, include_chain_arg = TRUE)
      
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
  
  # Return BOTH lists
  return(list(
    all_commands = all_commands_list,
    gpu_commands = gpu_command_list
  ))
}


#' Generate the core Python command string
.generate_python_command <- function(input_path_server, output_path_server, params, chain_index = NULL, include_chain_arg = TRUE) {
  transient_value <- ceiling(params$transient_proportion * params$samples * params$thin)
  verbose_value <- ceiling(params$samples / 10)
  
  python_args <- paste(
    "-m hmsc.run_gibbs_sampler", 
    "--input", shQuote(input_path_server), 
    "--output", shQuote(output_path_server),
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


#' Write the .sh scripts for Mode 3
.write_gpu_scripts <- function(gpu_command_list, run_config, output_script_dir, base_model_name) {
  for (gpu_index in 0:(run_config$gpu$n_gpus_available - 1)) {
    sh_filename <- file.path(output_script_dir, paste0(base_model_name, "_run_gpu_", gpu_index, ".sh"))
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
    
    con <- file(sh_filename, open = "wb")
    cat(c(header, wrapped_commands, footer), file = con, sep = "\n")
    close(con)
    
    message(paste("Successfully created sh script:", sh_filename))
  }
}

#' Takes the commands lists and writes the final .txt or .sh files.
#'
#' @param execution_mode The numeric mode (1, 2, or 3).
#' @param txt_commands A single character vector of all commands (for Modes 1/2).
#' @param gpu_commands A list of data frames (one per GPU) (for Mode 3).
#' @param output_script_dir Local path to save the script files.
#' @param base_model_name The *overall* model ID (e.g., "m_sp_tr_co_all")
#' @param run_config The full config list.
write_commands_scripts <- function(execution_mode, txt_commands, gpu_commands, 
                                     output_script_dir, base_model_name, run_config) {
  
  dir.create(output_script_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (execution_mode == 1 || execution_mode == 2) {
    # --- Write a single .txt log file for Modes 1 & 2 ---
    script_path <- file.path(output_script_dir, paste0(base_model_name, "_python_commands", ".txt"))
    message("\nGenerating master command log file: ", script_path)
    
    commands_to_write <- c(
      txt_commands,
      "",
      "# For manual run, remember to set the GPU:",
      "# CUDA_VISIBLE_DEVICES=0"
    )
    
    con <- file(script_path, open = "wb")
    cat(commands_to_write, file = con, sep = "\n")
    close(con)
    
    message(paste("Successfully created:", basename(script_path)))
    
  } else if (execution_mode == 3) {
    # --- Write separate .sh scripts for Mode 3 ---
    message("\nGenerating master .sh scripts for: ", base_model_name)
    
    .write_gpu_scripts(
      gpu_command_list = gpu_commands,
      run_config = run_config, 
      output_script_dir = output_script_dir, 
      base_model_name = base_model_name 
    )
  }
}

#' Create an Unfitted Hmsc Model for a cross-validation validation scheme
#'
#' Creates a new, unfitted `Hmsc` model object based on a training subset
#' of the data. This function is a core component of the k-fold 
#' cross-validation (CV) or spatial hold out (ho) workflow.
#'
#'
#' @param k [integer] The index (e.g., 1) of the fold to be *excluded* from this training model.
#' @param hM [Hmsc] The fitted, full `Hmsc` model object which will be used as a template.
#' @param partition [vector] A vector (same length as `hM$ny`) mapping each sampling unit 
#' to a fold (e.g., c(1, 2, 1, 1, 2, ...)). This is generated by `createPartition()`.
#'
#' @return [Hmsc] A new, *unfitted* `Hmsc` model object configured for
#'   the k-th training set.
#'
set_training_model <- function(k, hM, partition) {
  
  train_rows <- k != partition
  
  dfPi <- as.data.frame(matrix(NA,sum(train_rows),hM$nr),
                       stringsAsFactors = TRUE)
  colnames(dfPi) <- hM$rLNames
  dfPi <- droplevels(hM$dfPi[train_rows,, drop=FALSE])
  
  # for(r in seq_len(hM$nr)){
  #   dfPi[,r] <- factor(hM$dfPi[train_rows,r])
  # }
  
  switch(class(hM$X)[1L],
         matrix = {
           X_train <- hM$X[train_rows,,drop=FALSE]
         },
         list = {
           X_train <- lapply(hM$X, function(a) a[train_rows,,drop=FALSE])
         }
  )
  if(hM$ncRRR>0){
    XRRRTrain <- hM$XRRR[train_rows,,drop=FALSE]
  } else {
    XRRRTrain <- NULL
  }

  hM1 = Hmsc(Y=hM$Y[train_rows,,drop=FALSE], Loff=hM$Loff[train_rows,,drop=FALSE],
             XRRR=XRRRTrain, ncRRR = hM$ncRRR, XSelect = hM$XSelect,
             distr=hM$distr, studyDesign=dfPi, C=hM$C, ranLevels=hM$rL,
             XData = hM$XData[train_rows,,drop=FALSE], XFormula = hM$XFormula, 
             TrData = hM$TrData, TrFormula = hM$TrFormula
             )
  hM1 = setPriors(hM1, V0=hM$V0, f0=hM$f0, mGamma=hM$mGamma, UGamma=hM$UGamma, aSigma=hM$aSigma, bSigma=hM$bSigma, rhopw=hM$rhowp,
                  nfMax=3, nfMin=2)
  hM1$YScalePar = hM$YScalePar
  hM1$YScaled = (hM1$Y - matrix(hM1$YScalePar[1,],hM1$ny,hM1$ns,byrow=TRUE)) / matrix(hM1$YScalePar[2,],hM1$ny,hM1$ns,byrow=TRUE)
  hM1$XInterceptInd = hM$XInterceptInd
  hM1$XScalePar = hM$XScalePar
  switch(class(hM$X)[1L],
         matrix = {
           hM1$XScaled = (hM1$X - matrix(hM1$XScalePar[1,],hM1$ny,hM1$ncNRRR,byrow=TRUE)) / matrix(hM1$XScalePar[2,],hM1$ny,hM1$ncNRRR,byrow=TRUE)
         },
         list = {
           hM1$XScaled = list()
           for (zz in seq_len(length(hM1$X))){
             hM1$XScaled[[zz]] = (hM1$X[[zz]] - matrix(hM1$XScalePar[1,],hM1$ny,hM1$ncNRRR,byrow=TRUE)) / matrix(hM1$XScalePar[2,],hM1$ny,hM1$ncNRRR,byrow=TRUE)
           }
         }
  )
  if(hM1$ncRRR>0){
    hM1$XRRRScalePar = hM$XRRRScalePar
    hM1$XRRRScaled = (hM1$XRRR - matrix(hM1$XRRRScalePar[1,],hM1$ny,hM1$ncORRR,byrow=TRUE)) / matrix(hM1$XRRRScalePar[2,],hM1$ny,hM1$ncORRR,byrow=TRUE)
  }
  hM1$TrInterceptInd = hM$TrInterceptInd
  hM1$TrScalePar = hM$TrScalePar
  hM1$TrScaled = (hM1$Tr - matrix(hM1$TrScalePar[1,],hM1$ns,hM1$nt,byrow=TRUE)) / matrix(hM1$TrScalePar[2,],hM1$ns,hM1$nt,byrow=TRUE)
  
  return(hM1)
}

#' Import, Assemble, and Save a Single Hmsc gpu chain from a model
#'
#' This function processes a single MCMC run. It finds the posterior chain
#' files, imports them, assembles the fitted model, and saves the final .RData file.
#'
#' @param mcmc_params [data.frame row] A single row from the `run_config$mcmc`
#'   data frame, containing (samples, thin, n_chains, etc.).
#' @param run_nm [character] The unique name for this run (e.g., "fbs_M007_thin_1000...").
#' @param fitted_dir [character] The full file path where the final
#'   fitted model should be saved (e.g., "models/run_nm/fitted_run_nm.RData").
#'
#' @return Returns `invisible(TRUE)` on success, or `invisible(FALSE)` if the
#'   run was skipped or failed.
#'

import_posterior <- function(mcmc = mcmc_i, config = run_config, run_nm = run_name, 
                             partition_number = NULL, label = NULL){
  
  if(!is.null(label)){
    if(is.null(partition_number)) 
      stop("Argument label can be used only if a partition is load at the same time")
  }
  
  # --- A. have these runs already been imported?

  if(is.null(partition_number)){
    fitted_dir <- file.path("models", run_nm, paste0("fitted_", run_nm, ".rds"))
  }else{
    fitted_dir <- file.path("models", run_nm, label, paste0("fitted_", run_nm, "_", label, "_", partition_number, ".rds"))
  }
  
  if (file.exists(fitted_dir)){
    warning("Output file already exists. Returning path.")  
    return(fitted_dir)
  }
  
  # --- B. Find and Verify Posterior Chain Files ---
  run_specific_dir <- file.path("models", run_nm)
  
  if(!is.null(partition_number)){
    run_specific_dir <- file.path("models", run_nm, label)
  }
  
  if (!dir.exists(run_specific_dir)) {
    warning(paste("Directory not found for", run_nm, ". Skipping."))
    return(invisible(FALSE))
  }
  
  current_params <- mcmc
  chains <- current_params$n_chains
  chain_indices <- 0:(chains - 1)
  
  if(is.null(partition_number)){
    expected_filenames <- paste0(run_nm, "_post_chain", sprintf("%.2d", chain_indices), "_file.rds")
  }else{
    expected_filenames <- paste0(run_nm, "_", label, "_", partition_number,"_post_chain", sprintf("%.2d", chain_indices), "_file.rds")
  }
  
  chain_paths <- file.path(run_specific_dir, expected_filenames)
  
  existing_chains <- chain_paths[file.exists(chain_paths)]
  n_existing_chains <- length(existing_chains)
  
  if (length(existing_chains) < chains) {
    message(paste0("Found ", length(existing_chains), " of ", chains, " expected chains for ", run_nm, ". Take in mind!"))
  }
  
  # --- C. Import and Assemble the Model ---
  tryCatch({
    
    # Load the original unfitted model structure
    if(is.null(partition_number)){
      unfitted_path <- file.path("models", paste0("unfitted_", config$model_id, ".RData"))
      load(unfitted_path)
      unfitted_model <- models[[1]]  
    }else{
      unfitted_path <- file.path("models", run_nm, label, paste0("unfitted_", run_nm, "_", label, "_", partition_number, ".rds"))
      unfitted_model <- readRDS(unfitted_path)
    }
    
    
    post_list <- list()
    for (chain_path in existing_chains) {
      post_list[[length(post_list) + 1]] <- from_json(readRDS(chain_path)[[1]])[[1]]
    }
    
    samples <- current_params$samples
    thin <- current_params$thin
    transient <- ceiling(current_params$transient_proportion * samples * thin)
    
    fitted_model <- importPosteriorFromHPC(
      m = unfitted_model,
      postList = post_list,
      nSamples = samples,
      thin = thin,
      transient = transient,
      alignPost = TRUE
    )
    
    for(cInd in 1:n_existing_chains){
      if(unfitted_model$nr==1){
        for(s in 1:samples){
          fitted_model$postList[[cInd]][[s]]$Alpha = list(fitted_model$postList[[cInd]][[s]]$Alpha[1,])
        }
      }
      if(is.matrix(fitted_model$postList[[cInd]][[1]]$Alpha)){
        for(s in 1:samples){
          x = fitted_model$postList[[cInd]][[s]]$Alpha
          fitted_model$postList[[cInd]][[s]]$Alpha = lapply(seq_len(nrow(x)), function(k) x[k,])
        }
      }
    }
    
    # --- D. Save the Final, Individual Fitted Model Object ---
    saveRDS(fitted_model, file = fitted_dir)
    message("Successfully assembled and saved fitted model to: ", fitted_dir)
    
    if(is.null(partition_number)){
      return(invisible(TRUE))
    }else{
      return(fitted_dir)
    }
    
    
  }, error = function(e) {
    warning(paste("An error occurred while importing", run_nm, ":", e$message))
    return(invisible(FALSE))
  })
  
}

c.Hmsc <- getS3method("c","Hmsc")