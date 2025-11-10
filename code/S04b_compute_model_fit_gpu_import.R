# --- 1. Load Libraries and Utils ---
library(Hmsc)
library(jsonify)
source("code/config_model.R")
source("code/_utilities_hmsc_gpu.R")
set.seed(110724)

# --- 2. Configuration and Setup ---
paths <- list(
  local_dir = getwd(),
  models_dir = file.path(getwd(), "models"),
  unfitted_models_file = file.path("models", paste0("unfitted_", run_config$model_id, ".RData"))
)
dir.create(paths$models_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Define MCMC Parameters ---
mcmc_params <- run_config$mcmc
partition.sp <- NULL
Yc <- NULL
expected <- TRUE

# 4. Main Loop: Iterate Over the MCMC configuration
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  ptm <- proc.time() 
  # i <- 1
  parts <- 1:run_config$cv$k
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  mcmc_i <- run_config$mcmc[i, ]
  
  # --- A. Import to unfitted cross-validated model and read fitted cross-validation ---
  fitted_cv_paths <- lapply(
    X = parts, 
    FUN = function(X){
      import_posterior(
        mcmc = mcmc_i,
        run_nm = run_name,
        config = run_config, 
        partition_number = X
      )
    }
  ) 
  
  fitted_cv_models <- fitted_cv_paths |> 
    unlist() |> 
    lapply(readRDS)
  
  fitted_full_model_path <- file.path("models", run_name, paste0("fitted_", run_name, ".rds"))
  hM <- readRDS(fitted_full_model_path)
  postN <- Reduce(sum, lapply(hM$postList, length))
  
  # --- B. Call partition ---
  partition <- file.path(paths$models_dir, run_name[i], "cv", "partition.rds") |> 
    readRDS()
  
  predArray <- array(NA, c(hM$ny, hM$ns, postN))
  
  for (p in parts) {
    p <- 1
    message("predictions for partition ", p)
    val <- partition == p
    LoffVal = hM$LoffVal[val,,drop=FALSE]
    m <- fitted_cv_models[[p]]
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
    rm(pred1,val,m,postList,XData,dfPi)
    cat("Memory clean complete\n")
    cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
  }
  
  preds <- computePredictedValues(hM)
  cat("Calculating MF\n")
  MF <- evaluateModelFit(hM, predY=preds)
  cat("Calculating MFCV\n")
  MFCV <- evaluateModelFit(hM, predY=predArray)
  cat("Calculating WAIC\n")
  WAIC <- computeWAIC(hM)
  computational.time <- proc.time() - ptm
  cat("Time taken:", computational.time[3],"s \n\n")
  cat(sprintf("max mermory usage\n\tNcels: %.2f MB\n\tNcels: %.2f MB\n\tTotal: %.2f MB\n", gc()[11], gc()[12], gc()[11] + gc()[12]))
  
  dir_fit <- file.path("models", run_name, "model_fit") 
  dir.create(dir_fit, showWarnings = F, recursive = T)
  save(MF,MFCV,WAIC,predArray, preds, file = file.path(dir_fit, paste0("model_fit1", run_name, ".rds")))
  
}

