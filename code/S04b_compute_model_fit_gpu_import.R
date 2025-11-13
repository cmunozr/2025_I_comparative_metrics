# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
library(here)
source(file.path("code","config_model.R"))
source(file.path("code", "_utilities_hmsc_gpu.R"))
set.seed(110724)

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# This is a Crossvalidation model fit or a spatial hold-out
do_spatial_holdout <- FALSE

# --- 2. Configuration and Setup ---
models_dir <-  file.path(here::here(), "models")

# --- 3. Define MCMC Parameters ---
partition.sp <- NULL
Yc <- NULL
expected <- TRUE

# --- 4. Main Loop: Iterate Over the MCMC configuration
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  ptm <- proc.time() 
  # i <- 1
  
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  mcmc_i <- run_config$mcmc[i, ]
  
  fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  hM <- readRDS(fitted_full_model_path)
  postN <- Reduce(sum, lapply(hM$postList, length))
  
  if(do_spatial_holdout){
    label <- "ho"
    nfolds <- 1
    parts <- 1
  }else{
    label <- "cv"
    nfolds <- run_config$cv$k
    parts <- 1:nfolds
  }  
  
  # --- A. Import to unfitted validation models the posterior samples after GPU running at server
  # then read them as a fitted model ---
  fitted_val_paths <- lapply(
    X = parts, 
    FUN = function(X){
      import_posterior(
        mcmc = mcmc_i,
        run_nm = run_name,
        config = run_config, 
        partition_number = X,
        label = label
      )
    }
  ) 
  
  fitted_val_models <- fitted_val_paths |> 
    unlist() |> 
    lapply(readRDS)
  
  # --- B. Call partition ---
  partition <- file.path(models_dir, run_name[i], label, "partition.rds") |> 
    readRDS()
  
  # --- c. Make predictions for each part
  predArray <- array(NA, c(hM$ny, hM$ns, postN))
  
  for (p in parts) {
    # p <- 1
    message("predictions for partition ", p)
    val <- partition == p
    LoffVal = hM$LoffVal[val,,drop=FALSE]
    m <- fitted_val_models[[p]]
    m <- alignPosterior(m)
    postList <- poolMcmcChains(m$postList, start=1, thin=1)
    dfPi <- droplevels(hM$dfPi[val,, drop=FALSE])
    XDataVal <- hM$XData[val,, drop=FALSE]
    
    pred1 <- if (is.null(partition.sp)) {
      predict(m, post = postList, XData = XDataVal, 
                    XRRR = NULL, #hM$XRRR[val,, drop=FALSE],
                    Yc = NULL, #Yc[val,, drop=FALSE], 
                    studyDesign = dfPi,
                    mcmcStep = mcmcStep, expected = expected)
    } else {
      getSpeciesFoldPrediction(hM, m, val, postList, dfPi,
                               partition.sp = partition.sp,
                               mcmcStep = mcmcStep,
                               expected = expected,
                               nParallel = nParallel,
                               useSocket = useSocket)
    }
    cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
    cat("High memory use section, writing predictions to array\n")
    predArray[val,,] <- simplify2array(pred1)
    cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
    cat("Cleaning up memory\n")
    rm(pred1,val,m,postList,XDataVal,dfPi)
    cat("Memory clean complete\n")
    cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
  }
  
  #--- D. Calculating metrics
  
  preds <- computePredictedValues(hM)
  cat("Calculating MF\n")
  MF <- evaluateModelFit(hM, predY=preds)
  cat("Calculating MFCV\n")
  MFEVAL <- evaluateModelFit(hM, predY=predArray)
  cat("Calculating WAIC\n")
  WAIC <- computeWAIC(hM)
  computational.time <- proc.time() - ptm
  cat("Time taken:", computational.time[3],"s \n\n")
  cat(sprintf("max mermory usage\n\tNcels: %.2f MB\n\tNcels: %.2f MB\n\tTotal: %.2f MB\n", gc()[11], gc()[12], gc()[11] + gc()[12]))
  
  dir_fit <- file.path("models", run_name, "model_fit") 
  dir.create(dir_fit, showWarnings = F, recursive = T)
  
  saveRDS(list(MF = MF, MFEVAL = MFEVAL, WAIC = WAIC), file = file.path(dir_fit, paste0("metrics_", run_name, "_", label, nfolds, ".rds")))
  saveRDS(list(preds = preds, predArray = predArray), file = file.path(dir_fit, paste0("preds_", run_name, "_", label, nfolds, ".rds")))
  
}

