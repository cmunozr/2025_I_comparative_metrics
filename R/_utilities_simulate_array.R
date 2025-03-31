library(dplyr)

n.spp <- 40
n.sites <- 1000
samples <- 1000
  
spp_names <- openxlsx::int2col(1:n.spp)

spp <- runif(n.spp, min = 0, max = 1.5)


r_lognormPoi <- function(x, y){
  # Generate lognormal data
  lognormal_data <- rlnorm(x, 0, y)
  
  # Generate Poisson data using lognormal data as lambda
  poisson_data <- sapply(lognormal_data, function(x) rpois(1, x))
}


mpred_test <- lapply(spp, FUN = function(X){r_lognormPoi(x = n.sites, y = X)})
mpred_test <- do.call("cbind", mpred_test)
colnames(mpred_test) <- spp_names

traits_test <- matrix(nrow = n.spp, ncol = 3) |> 
  as.data.frame() |> 
  rename("wing.median" = V1, "kipps.median" = V2, "mass" = V3)

traits_test$wing.median <- runif(n.spp, min = 1, max = 6)
traits_test$kipps.median <- runif(n.spp, min = 1, max = 6)
traits_test$mass <- runif(n.spp, min = 1, max = 6)
rownames(traits_test) <- spp_names

mpre_test_list <- rep(list(mpred_test), samples)
# 1000 capas de 586 columnas y 10824 filas, Error: cannot allocate vector of size 23.6 Gb

filas <- nrow(mpre_test_list[[1]])
columnas <- ncol(mpre_test_list[[1]])
n_matrices <- length(mpre_test_list)

# Crear el array tridimensional vacío
arraytest <- array(NA, dim = c(filas, columnas, n_matrices),
                  dimnames = list(rownames(mpre_test_list[[1]]),
                                  colnames(mpre_test_list[[1]]),
                                  paste0(1:n_matrices)))

for (i in 1:n_matrices) {
  arraytest[,,i] <- mpre_test_list[[i]]
}

logic_arraytest <- arraytest > 0
storage.mode(logic_arraytest ) <- "logical"

binary_arraytest <- arraytest
storage.mode(binary_arraytest ) <- "integer"

rm("mpre_test_list", "mpred_test", "columnas", "filas", "n_matrices", "spp", "spp_names", "i")
gc()

