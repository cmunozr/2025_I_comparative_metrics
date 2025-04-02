# 10824 X 20 X 4000

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