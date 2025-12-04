# Setup
prediction_files <- list.files("results/predictions/metso", pattern = ".rds", full.names = TRUE)
output_dir <- "results/metrics/alpha_richness"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Loop through your batch files
for (f in prediction_files) {
  
  # 1. Read the batch (List of 4000)
  batch_list <- readRDS(f)
  
  # 2. Convert to Array [100, 67, 4000]
  batch_array <- simplify2array(batch_list)
  
  # 3. Calculate Metric: Alpha Richness
  # We want to sum across species (dimension 2) for every sample (dimension 3)
  # Result: Matrix [100 sites x 4000 samples]
  richness_matrix <- apply(batch_array, c(1, 3), sum)
  
  # 4. Save ONLY the metric
  # Extract batch number from filename for saving
  batch_name <- basename(f)
  saveRDS(richness_matrix, file = file.path(output_dir, batch_name))
  
  # Cleanup to free RAM
  rm(batch_list, batch_array, richness_matrix)
  gc()
}

# Richness
tm <- Sys.time()
a <- species_richness(pred.object = logic_predY)
tm_richness <- Sys.time()-tm
# 14.18454 secs

# geom abundance
tm <- Sys.time()
b <- geom_mean_abun(predY, a)
tm_geom<- Sys.time()-tm
# Time difference of 1.427866 mins

## Pre processing
traits_processed <- dbFD_preprocess_traits(x = TrData, a = logic_predY[,,1])
# CAMBIA SI USAMOS OTRA MATRIZ?

# Functional richness

tm <- Sys.time()
fr <- functional_richness(pred.object = logic_predY, trait.processed.object = traits_processed, stand.FRic = F, parallel = T)
tm_fr <- Sys.time()-tm


##	Community-weighted mean (CWM)

tm <- Sys.time()
CWM <- functcomp_parallel(pred.object = binary.predY, trait.processed.object = traits_processed, cwm.type = "dom", parallel = TRUE)
tm_CWM <- Sys.time()-tm

# RaoQ

tm <- Sys.time()
RaoQ <- divc_parallel(array_data = logic_predY, trait.processed.object = traits_processed, scale = FALSE, parallel = TRUE)
tm_RaoQ <- Sys.time()-tm

#--------------------------------------

# Richness
tm <- Sys.time()
a <- species_richness(pred.object = logic_arraytest)
a2 <- species_richness(pred.object = logic_arraytest2)
tm_richness2 <- Sys.time()-tm
    
# geom abundance
tm <- Sys.time()
b <- geom_mean_abun(arraytest, a)
tm_geom2<- Sys.time()-tm
    
# Pre processing
tm <- Sys.time()
traits_processed2 <- dbFD_preprocess_traits(x = traits_test, a = logic_arraytest[,,1])
tm_traits2<- Sys.time()-tm
    
# CAMBIA SI USAMOS OTRA MATRIZ?
  
# Functional richness
tm <- Sys.time()
fr <- functional_richness(pred.object = logic_arraytest, trait.processed.object = traits_processed2, stand.FRic = F, parallel = T)
tm_fr2 <- Sys.time()-tm
    
#	Community-weighted mean (CWM)
tm <- Sys.time()
CWM <- functcomp_parallel(pred.object = binary_arraytest[,,1:60], trait.processed.object = traits_processed2, cwm.type = "dom", parallel = T)
tm_CWM2 <- Sys.time()-tm
    
# RaoQ
tm <- Sys.time()
RaoQ <- divc_parallel(array_data = logic_arraytest, trait.processed.object = traits_processed2, scale = FALSE, parallel = TRUE)
tm_RaoQ2 <- Sys.time()-tm  

# no groups
tm <- Sys.time()
betas <- beta_diversity(pred.object = arraytest, group.each.rows = NULL)
tm_betas <- Sys.time()-tm  

# groups
    
# 10824 X 586 X 100
# tm_richness2: Time difference of 27.49811 secs
# tm_geom2: Time difference of 32.85858 secs
# tm_fr2: Time difference of 11.22873 mins
# tm_CWM2: Time difference of 10.45529 mins
# tm_RaoQ2: Time difference of 22.97023 mins
# tm_betas: 