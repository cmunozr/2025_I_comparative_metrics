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

# --- 3. Main Loop: Iterate Over the MCMC configuration
message(paste0("--- Import for Model ID: ", run_config$model_id))

for (i in 1:nrow(run_config$mcmc)) {
  # i <- 1
  run_name <- generate_run_name(run_config)[i]
  message(paste0("\n--- Processing Grid Row ", i, ": ", run_name, " ---"))
  
  if(do_spatial_holdout){
    label <- "ho"
    nfolds <- 1
    parts <- 1
  }else{
    label <- "cv"
    nfolds <- run_config$cv$k
    parts <- 1:nfolds
  }
  
  
  mf_object <- readRDS(file.path("models", run_name, "model_fit", paste0("metrics_", run_name, "_", label, nfolds, ".rds")))
  
  pdf(file = file.path("models", paste0(run_name, "_", label, nfolds, "_model_fit.pdf")))
  
  cMF = mf_object$MF
  cMFCV = mf_object$MFEVAL
  if(!is.null(cMF$TjurR2)){
    plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,
                     ": Tjur R2.\n",
                     "mean(MF) = ",as.character(mean(cMF$TjurR2,na.rm=TRUE)),
                     ", mean(MFCV) = ",as.character(mean(cMFCV$TjurR2,na.rm=TRUE))))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$R2)){
    plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,
                     ": R2. \n",
                     "mean(MF) = ",as.character(mean(cMF$R2,na.rm=TRUE)),
                     ", mean(MFCV) = ",as.character(mean(cMFCV$R2,na.rm=TRUE))))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$AUC)){
    plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,
                     ": AUC. \n",
                     "mean(MF) = ",as.character(mean(cMF$AUC,na.rm=TRUE)),
                     ", mean(MFCV) = ",as.character(mean(cMFCV$AUC,na.rm=TRUE))))
    abline(0,1)
    abline(v=0.5)
    abline(h=0.5)
  }
  if(FALSE && !is.null(cMF$O.TjurR2)){
    plot(cMF$O.TjurR2,cMFCV$O.TjurR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,": O.Tjur R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(FALSE && !is.null(cMF$O.AUC)){
    plot(cMF$O.AUC,cMFCV$O.AUC,xlim=c(0,1),ylim=c(0,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name, ": O.AUC"))
    abline(0,1)
    abline(v=0.5)
    abline(h=0.5)
  }      
  if(!is.null(cMF$SR2)){
    plot(cMF$SR2,cMFCV$SR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,
                     ": SR2. \n",
                     "mean(MF) = ",as.character(mean(cMF$SR2,na.rm=TRUE)),
                     ", mean(MFCV) = ",as.character(mean(cMFCV$SR2,na.rm=TRUE))))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }    
  if(FALSE && !is.null(cMF$C.SR2)){
    plot(cMF$C.SR2,cMFCV$C.SR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(run_name,": C.SR2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }  
  
  dev.off()
  
}
