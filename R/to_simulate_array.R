
spp_names <- openxlsx::int2col(1:586)

spp <- rpois(586, 10)

mpred_test <- lapply(spp, FUN = function(X){rpois(n = 10824, lambda = X)})
mpred_test <- do.call("cbind", mpred_test)
colnames(mpred_test) <- spp_names

traits_test <- matrix(nrow = 586, ncol = 3) |> 
  as.data.frame() |> 
  rename("wing.median" = V1, "kipps.median" = V2, "mass" = V3)

traits_test$wing.median <- runif(586, min = 1, max = 6)
traits_test$kipps.median <- runif(586, min = 1, max = 6)
traits_test$mass <- runif(586, min = 1, max = 6)
rownames(traits_test) <- spp_names

mpre_test_list <- rep(list(mpred_test), 100)
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

logic_arraytest <- arraytest
storage.mode(logic_arraytest ) <- "integer"

