#   This script performs the (preliminary) final post-processing of Joint Species Distribution 
#   Model (Hmsc) predictions. It calculates a suite of biodiversity metrics to 
#   compare two forestry scenarios: 
#     1. Conservation (METSO) 
#     2. Business-as-Usual (Control)
#
#   The script implements a "Cache and Load" pattern. It checks if a metric 
#   has already been calculated and saved to disk. If yes, it loads the object; 
#   if no, it performs the calculation and saves the result. This allows for 
#   efficient re-runs on local machines or HPC clusters.
#
# METRICS CALCULATED:
#   1. Species Richness (Alpha)
#   2. Mean Species Abundance (MSA) - [Schipper et al., 2020]
#   3. Potentially Disappeared Fraction (PDF) - [De Schryver et al., 2010]
#   4. Geometric Mean Abundance
#   5. Functional Richness (FRic)
#   6. Community-Weighted Mean (CWM)
#   7. Rao's Quadratic Entropy (RaoQ)
#   8. Sorensen Similarity

library(abind)
library(tidyverse)
library(Hmsc)
library(parallel)
library(here)

os <- Sys.info()['sysname']

if (os == "Windows") {
  cores <- 2
} else if (os %in% c("Linux", "Darwin")) {
  cores <- 10
}

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

# Setup

source(file.path("code", "diversity_metrics_functions.R"))
source(file.path("code", "config_model.R"))
output_dir <- file.path("results", "metrics")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
set.seed(11072024)

mfeval <- readRDS("models/fbs_M008_thin_1000_samples_1000_chains_4/model_fit/mfeval_fbs_M008_thin_1000_samples_1000_chains_4_ho1.rds")[["RMSE"]]
load("models/unfitted_fbs_M008.RData")
names(mfeval) <- colnames(models$fbs_M008$Y)

# WATCH, reduce array to work better to test
test <- T
expected_val <- T
if(expected_val){
  expected_string <- "_expected_true"
}else{
  expected_string <- "_expected_false"
}

# 1. Call data

if(test){
  num_spp <- 67
  samples <- 100
  to_sample <- list("a" = sample(seq(67), size = num_spp), "b"= 1:100)
}

pred_id_metso <- readRDS(file.path("results", paste0("predictions", expected_string), "metso", "pred_ids.rds"))
pred_id_control <- readRDS(file.path("results", paste0("predictions", expected_string), "control", "pred_ids.rds"))

matches <- readRDS(file.path("data", "metso", "raw", "matched_pairs.rds")) |> 
  mutate(in_pred_metso = standid_treated %in% pred_id_metso,
         in_pred_control = standid_matched_control %in% pred_id_control) |> 
  filter(in_pred_metso == 1)
  
# trait data
TrData <- readRDS(file.path("data", "traits", "TrData_complete.rds"))

# 2. Call predictions

sufixes <- c("metso", "control")

for(i in 1:length(sufixes)){
  # i <- 2
  sufix <- sufixes[i]
  prediction_files <- list.files(file.path("results", paste0("predictions", expected_string), sufix), pattern = "predY", full.names = TRUE)
  
  predY <- lapply(X = prediction_files, 
                  FUN = function(X){readRDS(X) |> simplify2array()})
  predY <- do.call("abind", c(predY, along = 1))
  
  if(sufix == "metso") predY <- predY[matches$in_pred_control, , ]
  
  if(i == 1) rownames(TrData) <- colnames(predY)
  
  ####>------- WATCH
  if(test) predY <- predY[ , to_sample[["a"]], to_sample[["b"]] ]
  ####>------- WATCH
  
  predY_logic <- predY > 0
  
  # dynamic name
  nm <- paste0("predY_", sufix)
  nm2 <- paste0("predY_logic_", sufix)
  
  assign(nm, predY)
  rm(predY); gc()
  assign(nm2, predY_logic)
  rm(predY_logic, sufix, prediction_files, nm, nm2); gc()
}

# 3. Calculate some metrics

cat(paste0("----- Writing to ", output_dir, "\n"))

# A. Richness
cat("----- Starting calculation of richness\n")
richness_file <- file.path("results", "metrics", "richness.rds")

if (!file.exists(richness_file)) {
  richness_metso <- species_richness(pred.object = predY_logic_metso)
  richness_control <- species_richness(pred.object = predY_logic_control)
  delta_richness <- richness_control - richness_metso 
  normalized <- scale(delta_richness)
  
  # Save all three so we can reload them later
  saveRDS(list(delta = delta_richness, norm = normalized, metso = richness_metso, control = richness_control), richness_file)
} else {
  cat("File exists, loading richness objects...\n")
  tmp <- readRDS(richness_file)
  # Load objects back into environment for dependent metrics (MSA, PDF, Geom)
  richness_metso <- tmp$metso
  richness_control <- tmp$control
  delta_richness <- tmp$delta
  rm(tmp)
}

# B. MSA
cat("----- Starting calculation of msa\n")
msa_file <- file.path("results", "metrics", "msa.rds")

if (!file.exists(msa_file)) {
  msa <- mean_species_abundance(pred.object.baseline = predY_metso, 
                                pred.object.scenario = predY_control, 
                                richness.baseline = richness_metso) 
  saveRDS(list(delta = msa, norm = scale(msa)), msa_file)
} else {
  cat("File exists, loading msa...\n")
  tmp <- readRDS(msa_file)
  msa <- tmp$delta
  rm(tmp)
}

# C. PDF
cat("----- Starting calculation of pdf\n")
pdf_file <- file.path("results", "metrics", "pdf.rds") 

if (!file.exists(pdf_file)) {
  pdf <- potentially_disappeared_fraction(richness.baseline = richness_metso, # Depends on A
                                          richness.scenario = richness_control,
                                          group.each.rows = 1)
  saveRDS(list(delta = pdf, norm = pdf), pdf_file)
} else {
  cat("File exists, loading pdf...\n")
  tmp <- readRDS(pdf_file)
  pdf <- tmp$delta
  rm(tmp)
}

# D. Geom abundance
cat("----- Starting calculation of geom\n")
geom_file <- file.path("results", "metrics", "geom.rds")

if (!file.exists(geom_file)) {
  geom_metso <- geom_mean_abun(pred.object = predY_metso, richness = richness_metso) # Depends on A
  geom_control <- geom_mean_abun(pred.object = predY_control, richness = richness_control) # Depends on A
  delta_geom <- geom_control - geom_metso 
  normalized <- scale(delta_geom)
  
  saveRDS(list(delta = delta_geom, norm = normalized, metso = geom_metso, control = geom_control), geom_file)
} else {
  cat("File exists, loading geom objects...\n")
  tmp <- readRDS(geom_file)
  geom_metso <- tmp$metso
  geom_control <- tmp$control
  delta_geom <- tmp$delta
  rm(tmp)
}

## Pre processing traits
# This is fast and required for the next steps, so we run it unconditionally
# to ensure TrData_processed exists for E, F, and G

fr_file <- file.path("results", "metrics", "fr.rds")
cwm_file <- file.path("results", "metrics", "CWM.rds")
raoq_file <- file.path("results", "metrics", "RaoQ.rds")

# LOGIC: Only calculate TrData_processed if we actually need to compute FR, CWM, or RaoQ
if (!file.exists(fr_file) || !file.exists(cwm_file) || !file.exists(raoq_file)) {
  
  cat("----- Starting calculation of TrData processing (Required for FR, CWM, or RaoQ)\n")
  
  # Define path for the processed traits cache
  tr_proc_file <- file.path("results", "metrics", "TrData_processed.rds")
  
  if (!file.exists(tr_proc_file)) {
    cat("----- Starting calculation of TrData processing (Heavy Step)\n")
    
    TrData_numeric <- TrData |> 
      select(where(is.numeric), -c("Total.individuals", "Female", "Male", "Unknown", "Complete.measures"))
    
    if(test) TrData_numeric <- TrData_numeric[to_sample[["a"]], ]
    
    # Calculate the heavy object
    TrData_processed <- dbFD_preprocess_traits(x = TrData_numeric, a = predY_logic_metso[,,1])
    
    # Save it immediately so we never have to run this again
    saveRDS(TrData_processed, tr_proc_file)
    
  } else {
    cat("----- Loading cached TrData_processed object...\n")
    TrData_processed <- readRDS(tr_proc_file)
  }
  
} else {
  cat("All trait-dependent metric files exist. Skipping TrData processing.\n")
}

# E. Functional richness
cat("----- Starting calculation of fr\n")
fr_file <- file.path("results", "metrics", "fr.rds")

if (!file.exists(fr_file)) {
  fr_metso <- functional_richness(pred.object = predY_logic_metso, trait.processed.object = TrData_processed, 
                                  stand.FRic = FALSE, parallel = TRUE, use.cores = cores)
  fr_control <- functional_richness(pred.object = predY_logic_control, trait.processed.object = TrData_processed, 
                                    stand.FRic = FALSE, parallel = TRUE, use.cores = cores)
  delta_fr <- fr_control - fr_metso
  normalized <- scale(delta_fr)
  
  saveRDS(list(delta = delta_fr, norm = normalized, metso = fr_metso, control = fr_control), fr_file)
} else {
  cat("File exists, loading fr objects...\n")
  tmp <- readRDS(fr_file)
  fr_metso <- tmp$metso
  fr_control <- tmp$control
  delta_fr <- tmp$delta
  rm(tmp)
}

# F. Community-weighted mean (CWM)
cat("----- Starting calculation of CWM\n")
cwm_file <- file.path("results", "metrics", "CWM.rds")

if (!file.exists(cwm_file)) {
  CWM_metso <- functcomp_parallel(pred.object = predY_logic_metso, trait.processed.object = TrData_processed, 
                                  cwm.type = "dom", parallel = TRUE, use.cores = 2)[[1]]
  CWM_control <- functcomp_parallel(pred.object = predY_logic_control, trait.processed.object = TrData_processed, 
                                    cwm.type = "dom", parallel = TRUE, use.cores = 2)[[1]]
  delta_CWM <- CWM_control - CWM_metso
  normalized <- scale(delta_CWM)
  
  saveRDS(list(delta = delta_CWM, norm = normalized, metso = CWM_metso, control = CWM_control), cwm_file)
} else {
  cat("File exists, loading CWM objects...\n")
  tmp <- readRDS(cwm_file)
  CWM_metso <- tmp$metso
  CWM_control <- tmp$control
  delta_CWM <- tmp$delta
  rm(tmp)
}

# G. RaoQ
cat("----- Starting calculation of RaoQ\n")
raoq_file <- file.path("results", "metrics", "RaoQ.rds")

if (!file.exists(raoq_file)) {
  RaoQ_metso <- divc_parallel(array_data = predY_logic_metso, trait.processed.object = TrData_processed, 
                              scale = FALSE, parallel = TRUE, use.cores = cores)
  RaoQ_control <- divc_parallel(array_data = predY_logic_control, trait.processed.object = TrData_processed, 
                                scale = FALSE, parallel = TRUE, use.cores = cores)
  
  # Ensure they are matrices before subtraction
  if(is.list(RaoQ_metso)) RaoQ_metso <- do.call(cbind, RaoQ_metso)
  if(is.list(RaoQ_control)) RaoQ_control <- do.call(cbind, RaoQ_control)
  
  delta_RaoQ <- RaoQ_control - RaoQ_metso
  normalized <- scale(delta_RaoQ)
  
  saveRDS(list(delta = delta_RaoQ, norm = normalized, metso = RaoQ_metso, control = RaoQ_control), raoq_file)
} else {
  cat("File exists, loading RaoQ objects...\n")
  tmp <- readRDS(raoq_file)
  RaoQ_metso <- tmp$metso
  RaoQ_control <- tmp$control
  delta_RaoQ <- tmp$delta
  rm(tmp)
}

# H. Sorensen similarity
cat("----- Starting calculation of sorensen\n")
sorensen_file <- file.path("results", "metrics", "sorensen.rds")

if (!file.exists(sorensen_file)) {
  sorensen <- sorensen_smilarity(pred.object.baseline = predY_metso, pred.object.scenario = predY_control)
  saveRDS(list(delta = sorensen, norm = scale(sorensen)), sorensen_file)
} else {
  cat("File exists, skipping sorensen calculation\n")
  tmp <- readRDS(sorensen_file)
  sorensen <- tmp$norm
  rm(tmp)
}

#I. Beta
# beta_metso <- beta_diversity(pred.object = predY_logic_metso, group.each.rows = 1)
# beta_control<- beta_diversity(pred.object = predY_logic_control, group.each.rows = 1)