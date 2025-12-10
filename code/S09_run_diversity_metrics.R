library(abind)
library(tidyverse)
library(Hmsc)

# Setup

source(file.path("code", "diversity_metrics_functions.R"))
source(file.path("code", "config_model.R"))
output_dir <- file.path("results", "metrics")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
set.seed(11072024)

####>------- WATCH, reduce array to work better to setup in laptop
num_spp <- 10
samples <- 100
to_sample <- list("a" = sample(seq(67), size = num_spp), "b"= sample(seq(4000), size = samples))
####>------- WATCH

# 1. Call data

pred_id_metso <- readRDS(file.path("results", "predictions", "metso", "pred_ids.rds"))
pred_id_control <- readRDS(file.path("results", "predictions", "control", "pred_ids.rds"))

matches <- readRDS(file.path("data", "metso", "donut_matches.rds")) |> 
  mutate(in_pred_metso = metso_standid %in% pred_id_metso,
         in_pred_control = control_standid %in% pred_id_control)
  
# trait data
TrData <- readRDS(file.path("data", "traits", "TrData_complete.rds"))

# 2. Call predictions

sufixes <- c("metso", "control")

for(i in 1:length(sufixes)){
  # i <- 2
  sufix <- sufixes[i]
  prediction_files <- list.files(file.path("results", "predictions", sufix), pattern = "predY", full.names = TRUE)
  
  predY <- lapply(X = prediction_files, 
                  FUN = function(X){readRDS(X) |> simplify2array()})
  predY <- do.call("abind", c(predY, along = 1))
  
  if(sufix == "metso") predY <- predY[matches$in_pred_control, , ]
  
  if(i == 1) rownames(TrData) <- colnames(predY)
  
  ####>------- WATCH
  predY <- predY[ , to_sample[["a"]], to_sample[["b"]] ]
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

# 3. Calculate

# A. Richness
richness_metso <- species_richness(pred.object = predY_logic_metso)
richness_control <- species_richness(pred.object = predY_logic_control)
delta_richness <- richness_metso - richness_control

# B. MSA
msa <- mean_species_abundance(pred.object.baseline = predY_metso, 
                              pred.object.scenario = predY_control, 
                              richness.baseline = richness_metso)

# C. PDF
pdf <- potentially_disappeared_fraction(richness.baseline = richness_metso, 
                                        richness.scenario = richness_control, 
                                        group.each.rows = 1)

# D. geom abundance
geom_metso <- geom_mean_abun(pred.object = predY_metso, richness = richness_metso)
geom_control <- geom_mean_abun(pred.object = predY_control, richness = richness_control)
delta_geom <- geom_metso - geom_control

## Pre processing traits

TrData_numeric <- TrData |> 
  select(where(is.numeric), -c("Total.individuals", "Female", "Male", "Unknown", "Complete.measures"))

####>------- WATCH
TrData_numeric <- TrData_numeric[to_sample[["a"]], ]
####>------- WATCH

TrData_processed <- dbFD_preprocess_traits(x = TrData_numeric, a = predY_logic_metso[,,1])

# E. Functional richness
fr_metso <- functional_richness(pred.object = predY_logic_metso, 
                                trait.processed.object = TrData_processed, 
                                stand.FRic = F, parallel = T)
fr_control <- functional_richness(pred.object = predY_logic_control, 
                                  trait.processed.object = TrData_processed, 
                                  stand.FRic = F, parallel = T)


##	Community-weighted mean (CWM)
tm <- Sys.time()
CWM <- functcomp_parallel(pred.object = binary.predY, trait.processed.object = traits_processed, cwm.type = "dom", parallel = TRUE)
tm_CWM <- Sys.time()-tm

# RaoQ
tm <- Sys.time()
RaoQ <- divc_parallel(array_data = logic_predY, trait.processed.object = traits_processed, scale = FALSE, parallel = TRUE)
tm_RaoQ <- Sys.time()-tm

