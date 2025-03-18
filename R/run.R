# 10824 X 20 X 4000

# Richness
tm <- Sys.time()
a <- species_richness()
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
tm_richness2 <- Sys.time()-tm

# geom abundance
tm <- Sys.time()
b <- geom_mean_abun(pred.object = arraytest, richness = a)
tm_geom2 <- Sys.time()-tm

traits_processed_test <- dbFD_preprocess_traits(x = traits_processed_test, a = arraytest, calc.FDiv = F, calc.FGR = F)
# pre procesar los rasgos toma 11 segundos en una potencial matriz de 10824 sitios por 586 especies

tm <- Sys.time()
fr <- functional_richness(pred.object = logic_arraytest, trait.object = traits_processed_test, stand.FRic = F, parallel = T)
tm_fr2 <- Sys.time() - tm

##	Community-weighted mean (CWM)

tm <- Sys.time()
CWM <- apply(X = arraytest, MARGIN = 3, FUN = function(X){FD::functcomp(x = traits_processed_test, a = X, CWM.type = "dom")})
tm_CWM2 <- Sys.time()-tm

# RaoQ

tm <- Sys.time()
RaoQ <- divc(array_data = logic_arraytest, trait.processed.object = traits_processed_test, scale = FALSE, parallel = TRUE)
tm_RaoQ <- Sys.time()-tm
