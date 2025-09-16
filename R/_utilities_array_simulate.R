library(dplyr)

set.seed(11072024)

n.spp <- 586
n.sites <- 1000
samples <- 100

spp_names <- openxlsx::int2col(1:n.spp)

spp <- runif(n.spp, min = 0, max = 1.5)

r_lognormPoi <- function(x, y) {
  lognormal_data <- rlnorm(x, 0, y)
  poisson_data <- sapply(lognormal_data, function(x) rpois(1, x))
  return(poisson_data)
}

# Crear el array tridimensional directamente
arraytest <- array(NA, dim = c(n.sites, n.spp, samples),
                   dimnames = list(1:n.sites, spp_names, paste0(1:samples)))

# Llenar el array con datos diferentes para cada muestra
for (i in 1:samples) {
  for (j in 1:n.spp) {
    arraytest[, j, i] <- r_lognormPoi(n.sites, spp[j])
  }
}
  
traits_test <- matrix(nrow = n.spp, ncol = 3) |> 
  as.data.frame() |> 
  rename("wing.median" = V1, "kipps.median" = V2, "mass" = V3)

traits_test$wing.median <- runif(n.spp, min = 1, max = 8)
traits_test$kipps.median <- runif(n.spp, min = 3, max = 4)
traits_test$mass <- rgamma(n.spp, shape = 3, scale = 1.5)
rownames(traits_test) <- spp_names

logic_arraytest <- arraytest > 0
storage.mode(logic_arraytest) <- "logical"

binary_arraytest <- logic_arraytest * 1
storage.mode(binary_arraytest) <- "integer"

rm("spp_names", "i", "n.spp", "samples", "n.sites", "r_lognormPoi", "spp", "j")
gc()

#-------
# reduction of abundance

reduce_abundance <- function(array.abundance, reduction.factor = 0.2, max.reduction = 8) {

  reduced_array <- array.abundance
  filas <- dim(array.abundance)[1]
  columnas <- dim(array.abundance)[2]
  muestras <- dim(array.abundance)[3]
  
  for (i in 1:filas) {
    for (j in 1:columnas) {
      for (k in 1:muestras) {
        abundance <- array.abundance[i, j, k]
        
        if (abundance > 4) {

          reduction_prob <- runif(1, min = 0, max = abundance * reduction.factor)
          
          if (reduction_prob > 1) {

            reduction <- sample(1:max.reduction, 1)
            
            reduced_array[i, j, k] <- max(0, abundance - reduction)
          }
        }
      }
    }
  }
  storage.mode(reduced_array) <- "integer"
  return(reduced_array)
}

arraytest2 <- reduce_abundance(arraytest, reduction.factor = 1, max.reduction = 15)

logic_arraytest2 <- arraytest2 > 0
storage.mode(logic_arraytest2) <- "logical"

binary_arraytest2 <- logic_arraytest2 * 1
storage.mode(binary_arraytest2) <- "integer"
