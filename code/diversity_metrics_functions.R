# arguments names with .
# function names with _
# Objects names _

# rm(list=setdiff(ls(), c("predY", "log_zero", "TrData", "logic_predY", "a", "b", "traits_processed_test", "logic_arraytest", "binary.predY", "traits_processed"))); gc()

# Metrics in management

##  MSA, mean species abundance
## 10.1111/gcb.14848
## "the calculation of the MSA metric, where IAR denotes individual species' abundance in an undisturbed 
## reference situation, IAI the abundance of the species in the impacted situation, and IAI/IAR the 
## truncated abundance ratio, which is calculated only for original species (i.e. occurring in the reference situation)""

mean_species_abundance <- function(pred.object.baseline, pred.object.scenario, richness.baseline, group.each.rows = NULL){
  
  if (!identical(dim(pred.object.baseline), dim(pred.object.scenario))) {
    stop("Dimensions of 'pred.object.baseline' and 'pred.object.scenario' must be identical.")
  }
  
  if (!identical(dim(pred.object.baseline)[c(1,3)], dim(richness.baseline))) {
    stop("Dimensions of 'richness.baseline' [sites, samples] do not match site/sample dimensions of prediction arrays.")
  }
  
  # ratio of abundance between posterior predictions (AI/AR)
  ratio <- (pred.object.scenario / pred.object.baseline)
  # as we are not predicting in the same place of training, there is no a 
  # 1:1 match between both sets
  
  # handle Inf: it means value of baseline is 0, so the species was not present in the baseline, according to definition of MSA should be 0
  ratio[is.infinite(ratio)] <- 0
  
  # and NaN: it means value of baselina and prediction is 0, so they shoudl be excluded
  ratio[is.nan(ratio)] <- NA
  
  # Truncated to 1
  ratio[ratio > 1] <- 1
  
    # Sum ratios of change on sites of each replicate of the posterior
  sum_ratios <- apply(ratio, c(1,3), sum, na.rm = T)
  
  # Divide the sum of ratios between the number of species per site (richness previously calculated)
  msa <- sum_ratios/richness.baseline
  
  # If richness.baseline was 0, result is NaN or Inf. 
  # Convert this to NA, as MSA is undefined for empty baselines.
  msa[is.nan(msa) | is.infinite(msa)] <- NA
  
  return(msa)
  
  ###### MISSING mean grouped
}

#------------------------

##  PDF, potentially disappeared fraction
## 10.1007/s11367-010-0205-2
## "the species number on the baseline land use type (Sb) divided by and the species number on the occupied land use type i (Si)"

potentially_disappeared_fraction <- function(richness.baseline, richness.scenario, group.each.rows = NULL) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  library(dplyr)
  
  if (!identical(dim(richness.baseline), dim(richness.scenario))) {
    stop("Dimensions of 'richness.baseline' and 'richness.scenario' must be identical.")
  }
  
  N_sites_total <- nrow(richness.baseline)
  N_posterior <- ncol(richness.baseline)
  
  # Determine number of groups based on rows (sites)
  if (!is.null(group.each.rows) && group.each.rows > 0 && group.each.rows < N_sites_total) {
    groups <- ceiling(N_sites_total / group.each.rows)
    group.each.rows <- as.integer(group.each.rows)
  } else {
    groups <- 1
    group.each.rows <- N_sites_total
  }
  
  # Calculate local CF

  ratio <- 1 - (richness.scenario/richness.baseline)
  
  ratio[ is.infinite(ratio) | is.nan(ratio) ] <- NA

  # to save
  results_list <- vector("list", groups)
  
  for (g in 1:groups) {
    # g <- 1
    # Calculate site indices for the current group
    r0 <- (g - 1) * group.each.rows + 1
    r1 <- min(g * group.each.rows, N_sites_total)
    current_site_indices <- r0:r1
    n_group_sites <- length(current_site_indices)
    
    ratio_group <- ratio[current_site_indices, , drop = FALSE]
    
    mean_group <- colMeans(ratio_group, na.rm = T) |> 
      t() |> 
      as.data.frame() # only possible as metso and control have similar area and when grouping dont affect it is setuped to 1###### WATCH
    
    results_list[[g]] <- mean_group
    
  } # End group loop
  
  m.res <- dplyr::bind_rows(results_list)
  
  return(m.res)
}

## STARTt threat abatement component of STAR
## https://doi.org/10.1038/s41559-021-01432-0


# STARti = sum((Ps * Ws, Cst)/Ni)
# Ps = current AOH of species s within location i (as a percentage of the global AOH)
# Ws = IUCN Red List category weight of species s
# Cst = relative contribution of threat t to the extinction risk of species s
# Ni = total number of species at location i



#library(rredlist)

#rl_species("Gorilla", "gorilla")

# From ecology 

#------------

##	Species richness

species_richness <- function(pred.object = logic_predY){
  
  dm <- dim(pred.object)
  
  m_res <- matrix(nrow = dm[1], ncol = dm[3])
  
  for(i in 1:dm[3]){
    m_res[ , i] <- pred.object[,,i] |> 
      rowSums()
  }
  
  storage.mode(m_res) <- "integer"
  
  return(m_res)
}

#-------------------

##	Geometric mean abundance
## https://doi.org/10.1016/j.biocon.2016.08.024

geom_mean_abun <- function(pred.object = predY, richness = a){
  
  dm <- dim(pred.object)
  
  m_res <- matrix(nrow = dm[1], ncol = dm[3])
  
  for(i in 1:dm[3]){
    tmp <- pred.object[,,i] |> 
      log_zero() |>
      rowSums( na.rm = T)/richness[,i]
    m_res[ , i] <- exp(tmp) |> 
      round(3)
  }
  return(m_res)
  
}

#-----------------

## Sorensen similarity
## https://doi.org/10.1016/j.biocon.2016.08.024

sorensen_smilarity <- function(pred.object.baseline, pred.object.scenario){
  
  if (!identical(dim(pred.object.baseline), dim(pred.object.scenario))) {
    stop("Dimensions of 'pred.object.baseline' and 'pred.object.scenario' must be identical.")
  }
  
  sum_abs_dif <- (pred.object.scenario - pred.object.baseline) |> 
    abs() |> 
    apply(3, rowSums)
  
  total_abundance <- (pred.object.scenario + pred.object.baseline) |> 
    abs() |> 
    apply(3, rowSums)
  
  # handle Inf and NaN
  relative_abundance <- sum_abs_dif/total_abundance
  
  relative_abundance <- ifelse(is.infinite(relative_abundance) | is.nan(relative_abundance), NA, relative_abundance)
  
  return(1 - relative_abundance)
  
}

#------------------

##	Functional richness 
##  https://doi.org/10.1890/08-2244.1

functional_richness <- function(pred.object, trait.processed.object, stand.FRic, parallel = T, use.cores = 2) {
  
  # time process array with 100 matrix of 10824 rows x 20 cols: Parallel = 7.814091 mins, No Parallel = 20.49533 mins
  # time process array with 100 matrix of 10824 rows x 586 cols: Parallel = 12.7466 mins
  
  if (!requireNamespace("geometry", quietly = TRUE)) install.packages("geometry")
  if (!requireNamespace("FD", quietly = TRUE)) install.packages("FD")
  
  library(geometry)
  library(parallel)
  library(FD)

  
  # Calculate functional diversity for a matrix dxm
  # - each site is a row
  # - columns are species
  
  process_matrix <- function(c) {
    
    FRic <- numeric(dm[1])
    
    for (i in seq_len(dm[1])){
      sppres <- pred.object[i, , c]
      S <- sum(sppres)
      nb.sp <- S
      tr.FRic <- data.frame(traits[sppres, ])
      if(S < 4){
        FRic[i] <- NA
        next
      }
      if (all(x.class2 == "factor" | x.class2 == "ordered")) {
        if (length(x.class2) == 1 & x.class2[1] == "ordered") {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
        else {
          if (!stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))
          if (stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))/FRic.all
        }
      }
      else {
        if (dim(tr.FRic)[2] > 1 & nb.sp >= 3) {
          if (war) 
            thresh <- 4
          if (!war) 
            thresh <- 3
          if (nb.sp >= thresh) {
            cvhull <- geometry::convhulln(tr.FRic, "FA")
            if (!stand.FRic) 
              FRic[i] <- cvhull$vol
            if (stand.FRic) 
              FRic[i] <- cvhull$vol/FRic.all
          }
          else {
          }
        }
        if (dim(tr.FRic)[2] == 1) {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
      }
    }
    return(round(FRic, 3))
  }
  
  traits <- trait.processed.object$traits.FRic
  x.class2 <- trait.processed.object$x.class2
  hull.all <- trait.processed.object$hull.all
  FRic.all <- trait.processed.object$FRic.all
  war <- trait.processed.object$warning
  
  species_names <- colnames(pred.object)
  species_index <- setNames(seq_along(species_names), species_names)
  
  colnames(pred.object) <- species_index[colnames(pred.object)]
  rownames(traits) <- as.character(species_index[rownames(traits)])
  trait_list <- split(traits, rownames(traits))
  
  dm <- dim(pred.object)
  
  m.res <- matrix(nrow = dm[1], ncol = dm[3])
  
  if(parallel){
    
    num_cores <- use.cores
    
    if(.Platform$OS.type == "windows"){
      
      cl <- makeCluster(num_cores)
      results <- parLapply(cl, seq_len(dm[3]), process_matrix)
      stopCluster(cl)
      
    } else {
      results <- mclapply(seq_len(dm[3]), process_matrix, mc.cores = num_cores)
    }
    
  } else {
    results <- lapply(seq_len(dm[3]), process_matrix)
  }
  
  # Combinar los resultados en la matriz m.res
  for (c in seq_len(dm[3])) {
    m.res[, c] <- results[[c]]
  }
  
  return(m.res)
}

#-----------------------------

##  Rao’s quadratic entropy coefficient (RaoQ)   
##  https://doi.org/10.1890/08-2244.1

divc_parallel <- function(array_data, trait.processed.object = NULL, scale = FALSE, parallel = FALSE, use.cores = 2) {
  
  if (!requireNamespace("FD", quietly = TRUE)) install.packages("FD")
  library(FD)
  
  div_calc <- function(matrix_data) {
    
    df <- as.data.frame(t(matrix_data))
    
    if (!inherits(df, "data.frame")) 
      stop("Non convenient df")
    
    if (any(df < 0)) {
      stop("Negative value in df")
    }
    
    if (!is.null(dis)) {
      if (!inherits(dis, "dist")) {
        stop("Object of class 'dist' expected for distance")
      }
      dis_matrix <- as.matrix(dis)
      if (nrow(df) != nrow(dis_matrix)) {
        stop("Non convenient df")
      }
      dis_dist <- as.dist(dis_matrix)
    } else {
      dis_dist <- as.dist((matrix(1, nrow(df), nrow(df)) - diag(rep(1, nrow(df)))) * sqrt(2))
    }
    
    results <- sapply(1:ncol(df), function(i) {
      if (sum(df[, i]) < 1e-16) {
        return(0)
      } else {
        return((t(df[, i]) %*% (as.matrix(dis_dist)^2) %*% df[, i]) / (2 * (sum(df[, i])^2)))
      }
    })
    return(results)
  }
  
  dis <- trait.processed.object$x.dist
  
  if (parallel) {
    num_cores <- use.cores
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(num_cores)
      #parallel::clusterExport(cl, c("array_data", "dis", "scale"))
      results <- parallel::parApply(cl, array_data, 3, div_calc)
      parallel::stopCluster(cl)
    } else {
      n_samples <- dim(array_data)[3]
      
      # Run mclapply over the indices, extracting the matrix slice [,,i] for each core
      results <- parallel::mclapply(1:n_samples, function(i) {
        div_calc(array_data[,,i]) 
      }, mc.cores = num_cores)
      results <- do.call(cbind, results)
    }
  } else {
    results <- apply(array_data, 3, div_calc)
  }
  
  # Escalar si es necesario
  if (scale) {
    divmax <- divcmax(dis)$value
    results <- results / divmax
  }
  
  return(results)
}

#----------------------

## CWM: community weighted mean
## https://doi.org/10.1890/08-2244.1

functcomp_parallel <- function(pred.object, trait.processed.object, cwm.type = "dom", parallel = FALSE, use.cores = 2) {
  
  require(FD)
  require(dplyr)
  require(parallel)
  
  traits <- trait.processed.object$traits.FRic
  num_traits <- ncol(traits)
  num_sites <- dim(pred.object)[1]
  num_pred <- dim(pred.object)[3]
  
  # Inicializar listas para almacenar las columnas de cada rasgo
  trait_lists <- vector("list", num_traits)
  
  process_matrix <- function(i) {
    trait_weightead <- FD::functcomp(x = trait_nt, a = pred.object[,,i], CWM.type = cwm.type)
    return(trait_weightead)
  }
  
  for(nt in 1:num_traits){
    
    trait_nt <- traits |> select(all_of(nt))
    
    if (parallel) {
      num_cores <- use.cores
      if (.Platform$OS.type == "windows") {
        cl <- makeCluster(num_cores)
        clusterExport(cl, c("trait_nt", "pred.object", "cwm.type"), envir = environment())
        sublist <- parLapply(cl, seq_len(num_pred), process_matrix)
        stopCluster(cl)
      } else {
        sublist <- mclapply(seq_len(num_pred), process_matrix, mc.cores = num_cores)
      }
    } else {
      sublist <- lapply(seq_len(num_pred), process_matrix)
    }
    
    
    trait_lists[[nt]] <- do.call("cbind", sublist)
    rm("sublist")
  }
  
    
  return(trait_lists)
}

#-----------------------------

##	Beta: https://onlinelibrary-wiley-com.ludwig.lub.lu.se/doi/full/10.1111/ele.12141
##  VarY the variation in species composition among sites in the geographical region of interest
## LCBD contributions of individual sampling units to the overall beta diversity
## SCBD contributions of individual species to the overall beta diversity
## Agrupamiento no espacial
## Agrupamiento espacial(?)
## Paralleling (?)

# value
# final_varY_df: sample, group, varY (BDTotal for the group)
# final_LCDB_df: sample, group, site_global_idx, LCDB (within the group)
# final_SCDB_df: sample, group, species_idx, SCDB (within the group)

beta_diversity <- function(pred.object = arraytest, group.each.rows = NULL){

  require(dplyr)
  
  dm <- dim(pred.object)
  N_samples <- dm[3]
  N_sites_total <- dm[1]
  N_species <- dm[2]
  
  # Determine number of groups
  if(!is.null(group.each.rows) && group.each.rows > 0 && group.each.rows < N_sites_total){
    groups <- ceiling(N_sites_total / group.each.rows)
  } else {
    groups <- 1
    group.each.rows <- N_sites_total # Ensure full range if groups=1
  }

  varY_list <- vector("list", N_samples * groups)
  LCDB_list <- vector("list", N_samples * groups)
  SCDB_list <- vector("list", N_samples * groups)
  
  list_idx <- 1 # Index for storing results
  
  
  for(i in 1:N_samples) { 
  # i <- 2
    for(g in 1:groups) {
      # g <- 1
      # Site indices for current group
      r0 <- (g - 1) * group.each.rows + 1
      r1 <- min(g * group.each.rows, N_sites_total)
      current_site_indices <- r0:r1
      n_group_sites <- length(current_site_indices)
      
      # Extract group submatrix
      Y_group <- pred.object[current_site_indices, , i, drop = FALSE] 
      
      # Calculations for the group
      rsums <- rowSums(Y_group)
      denom <- ifelse(rsums == 0, 1, rsums) # Avoid division by zero
      hellinger_transf <- sqrt(Y_group / denom) 
      
      # Explicitly index rows to set to zero
      rows_to_zero <- which(rsums == 0)
      if (length(rows_to_zero) > 0) {
        print(g)
        hellinger_transf[rows_to_zero, ] <- 0 
      }
      
      # Calculate based on group's centroid
      col_means_H_group <- colMeans(hellinger_transf) 
      squared_deviations_site <- sweep(hellinger_transf, MARGIN = 2, STATS = col_means_H_group, FUN = "-")^2
      SS <- sum(squared_deviations_site) # Group's SSTotal
      
      varY_group <- if (n_group_sites > 1) SS / (n_group_sites - 1) else NA
      
      # Group's LCBD & SCBD - handle SS=0 case
      if (SS > 1e-12) { # Use tolerance for floating point
        LCDB_group <- rowSums(squared_deviations_site) / SS 
        SCDB_group <- colSums(squared_deviations_site) / SS 
      } else {
        LCDB_group <- rep(NA, n_group_sites)
        SCDB_group <- rep(NA, N_species)
      }
      
      varY_list[[list_idx]] <- tibble(sample = i, group = g, varY = varY_group)
      LCDB_list[[list_idx]] <- tibble(sample = i, group = g, site_global_idx = current_site_indices, LCDB = LCDB_group)
      SCDB_list[[list_idx]] <- tibble(sample = i, group = g, species_idx = 1:N_species, SCDB = SCDB_group)
      
      list_idx <- list_idx + 1 # Increment list index
      
    } 
    
  }
  
  varY <- dplyr::bind_rows(varY_list) |> as.data.frame()
  LCDB <- dplyr::bind_rows(LCDB_list) |> as.data.frame()
  SCDB <- dplyr::bind_rows(SCDB_list) |> as.data.frame()
  
  return(list("varY" = varY, "LCDB" = LCDB, "SCDB" = SCDB))
  
}

#-----------------------------

# utilities

log_zero <- function(x) {
  ifelse(x == 0, 0, log(x))
}

# preprocess traits. Taken from FD function dbFD (original function quite slowly)
##  https://doi.org/10.1890/08-2244.1

dbFD_preprocess_traits <- function (x, a, w, w.abun = TRUE, stand.x = TRUE, ord = c("podani", "metric"), 
                                    asym.bin = NULL, corr = "sqrt", calc.FRic = TRUE, m = "max", stand.FRic = FALSE, 
                                    scale.RaoQ = FALSE, calc.FGR = F, clust.type = "ward", 
                                    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, 
                                    km.crit = c("calinski", "ssi"), calc.CWM = TRUE, CWM.type = c("dom", "all"), 
                                    calc.FDiv = F, dist.bin = 2, print.pco = FALSE, messages = TRUE) 
{
  # w.abun = TRUE; stand.x = TRUE; ord = c("podani", "metric"); asym.bin = NULL;
  # corr = c("sqrt", "cailliez", "lingoes", "none"); calc.FRic = TRUE; m = "max"; 
  # stand.FRic = FALSE; scale.RaoQ = FALSE; calc.FGR = FALSE; clust.type = "ward"; 
  # km.inf.gr = 2; km.sup.gr = nrow(x) - 1; km.iter = 100; km.crit = c("calinski", "ssi"); 
  # calc.CWM = TRUE; CWM.type = c("dom", "all"); calc.FDiv = TRUE; dist.bin = 2; 
  # print.pco = FALSE; messages = TRUE
  
  if (!requireNamespace("geometry", quietly = TRUE)) install.packages("geometry")
  if (!requireNamespace("FD", quietly = TRUE)) install.packages("FD")
  
  library(geometry)
  library(parallel)
  library(FD)
  
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  CWM.type <- match.arg(CWM.type)
  km.crit <- match.arg(km.crit)
  
  tol <- .Machine$double.eps
  if (!is.logical(messages)) 
    stop("'messages' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.FRic)) 
    stop("'stand.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.x)) 
    stop("'stand.x' must be TRUE or FALSE.", "\n")
  if (!is.logical(w.abun)) 
    stop("'w.abun' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FRic)) 
    stop("'calc.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FDiv)) 
    stop("'calc.FDiv' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FGR)) 
    stop("'calc.FGR' musts be TRUE or FALSE.", "\n")
  if (!is.logical(calc.CWM)) 
    stop("'calc.CWM' must be TRUE or FALSE.", "\n")
  if (!is.logical(scale.RaoQ)) 
    stop("'scale.RaoQ' must be TRUE or FALSE.", "\n")
  if (!is.logical(print.pco)) 
    stop("'print.pco' must be TRUE or FALSE.", "\n")
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- dim(x)[1]
    t.x <- dim(x)[2]
    if (is.null(row.names(x))) 
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x))) 
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels"))) 
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (missing(a)) {
    ab.names <- list("Community1", x.rn)
    a <- matrix(1, 1, s.x, dimnames = ab.names)
  }  else {
    if (is.matrix(a) | is.data.frame(a)) {
      s.a <- dim(a)[2]
      ab.t <- t(a)
      if (is.null(row.names(ab.t))) 
        stop("'a' must have column names.", "\n")
      else ab.t.row <- row.names(ab.t)
      a <- as.matrix(a)
    }
    if (is.vector(a)) {
      s.a <- length(a)
      if (is.null(names(a))) 
        stop("'a' must have names.", "\n")
      else ab.t.row <- names(a)
      ab.names <- list("Community1", ab.t.row)
      a <- matrix(a, 1, s.a, dimnames = ab.names)
    }
    if (s.x != s.a) 
      stop("Different number of species in 'x' and 'a'.", 
           "\n")
    if (any(ab.t.row != x.rn)) 
      stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
           "\n")
  }
  
  a <- as.matrix(a)
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    stop("At least one community has zero-sum abundances (no species).", 
         "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", 
         "\n")
  if (!missing(w) & is.dist.x) 
    stop("When 'x' is a distance matrix, 'w' should be left missing.", 
         "\n")
  if (!missing(w) & !is.dist.x) {
    if (!is.numeric(w) | length(w) != t.x) 
      stop("'w' should be a numeric vector of length = number of traits.", 
           "\n")
    else w <- w/sum(w)
  }
  if (missing(w)) 
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) 
        x[, x.class == "character"] <- as.factor(x[, 
                                                   x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1])) 
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          x.rn <- x.rn[-pos.NA]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)])) 
            stop("'dist.bin' must be an integer between 1 and 10.", 
                 "\n")
          x.dist <- dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = T, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.", 
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x))) 
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).", 
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist))) 
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0)) 
      stop("At least one species has no trait data.", "\n")
  }
  c <- dim(a)[1]
  if (!w.abun) 
    for (h in 1:c) {
      abpos <- which(a[h, ] > 0)
      a[h, abpos] <- 1
    }
  attr(x.dist, "Labels") <- x.rn
  if (is.euclid(x.dist)) 
    x.dist2 <- x.dist
  if (!is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- lingoes(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", 
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", 
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2)) 
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", 
             "\n")
      if (is.euclid(x.dist2)) 
        if (messages) 
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", 
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages) 
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", 
            "\n")
    }
  }
  x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  traits <- round(x.pco$li, .Machine$double.exponent)
  nb.sp <- numeric(c)
  for (i in 1:c) {
    sp.pres <- which(a[i, ] > 0)
    traits.sp.pres <- traits[sp.pres, , drop = F]
    traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < 
                     tol] <- 0
    nb.sp[i] <- nrow(unique(traits.sp.pres))
  }
  names(nb.sp) <- row.names(a)
  min.nb.sp <- min(nb.sp)
  if (min.nb.sp < 3) 
    if (messages) 
      cat("FEVe: Could not be calculated for communities with <3 functionally singular species.", 
          "\n")
  if (min.nb.sp < 2) 
    if (messages) 
      cat("FDis: Equals 0 in communities with only one functionally singular species.", 
          "\n")
  if (calc.FRic) {
    x.class2 <- sapply(x, data.class)
    if (all(x.class2 == "factor" | x.class2 == "ordered")) {
      if (length(x.class2) == 1 & x.class2[1] == "ordered") {
        traits.FRic1 <- rank(x[, 1])
        names(traits.FRic1) <- x.rn
        traits.FRic <- data.frame(traits.FRic1)
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one ordinal trait present in 'x'. FRic was measured as the range of the ranks, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when 'x' is a single ordinal trait.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      else {
        traits.FRic <- x
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only categorical and/or ordinal trait(s) present in 'x'. FRic was measured as the number of unique trait combinations, NOT as the convex hull volume.", 
              "\n")
        if (stand.FRic) 
          FRic.all <- nrow((unique(traits.FRic)))
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when only categorical and/or ordinal trait(s) present in 'x'.", 
                "\n")
        }
      }
    }
    else {
      if (x.pco$nf == 1) {
        traits.FRic <- x.pco$li
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      if (x.pco$nf > 1) {
        warning <- FALSE
        m.max <- min.nb.sp - 1
        if (m == "min") {
          warning <- TRUE
          if (min.nb.sp < 4) {
            nb.sp2 <- nb.sp[nb.sp > 3]
            m.min <- floor(log2(min(nb.sp2)))
            if (messages) 
              cat("FRic: To respect s >= 2^t, FRic could not be calculated for communities with <4 functionally singular species.", 
                  "\n")
          }
          else m.min <- floor(log2(min.nb.sp))
        }
        else {
          if (min.nb.sp < 3) {
            nb.sp2 <- nb.sp[nb.sp > 2]
            m.max <- min(nb.sp2) - 1
            if (messages) 
              cat("FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species.", 
                  "\n")
          }
          else m.max <- m.max
        }
        if (is.numeric(m) & m <= 1) 
          stop("When 'm' is an integer, it must be >1.", 
               "\n")
        if (is.numeric(m) & m > m.max) 
          m <- m.max
        if (m == "min") 
          m <- m.min
        if (m == "max") 
          m <- m.max
        if (!is.numeric(m) & m != "min" & m != "max") 
          stop("'m' must be an integer >1, 'min', or 'max'.", 
               "\n")
        if (m < x.pco$nf) {
          traits.FRic <- x.pco$li[, 1:m]
          if (x.pco$nf - m == 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last PCoA axis (out of", 
                  x.pco$nf, "in total) was removed.", "\n")
          if (x.pco$nf - m > 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last", 
                  x.pco$nf - m, "PCoA axes (out of", x.pco$nf, 
                  "in total) were removed.", "\n")
          if (is.euclid(x.dist)) {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr != "none") {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (based on corrected distance matrix) =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr == "none") {
            delta <- -0.5 * bicenter.wt(x.dist * x.dist)
            lambda <- eigen(delta, symmetric = TRUE, 
                            only.values = TRUE)$values
            sum.m <- sum(lambda[1:m])
            sum.n <- sum(lambda)
            lambda.neg <- c(lambda[lambda < 0])
            max.neg <- abs(min(lambda.neg))
            qual.FRic <- (sum.m + (length(lambda[1:m]) * 
                                     max.neg))/(sum.n + ((length(lambda) - 1) * 
                                                           max.neg))
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) =", 
                  qual.FRic, "\n")
          }
        }
        if (m >= x.pco$nf) {
          qual.FRic = 1
          traits.FRic <- x.pco$li
          if (x.pco$nf == 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. The 2 PCoA axes were kept as 'traits'.", 
                  "\n")
          if (x.pco$nf > 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. All", 
                  x.pco$nf, "PCoA axes were kept as 'traits'.", 
                  "\n")
        }
        hull.all <- convhulln(traits.FRic, "FA")
        FRic.all <- hull.all$vol
      }
    }
  }
  res <- list()
  res$traits.FRic <- traits.FRic
  res$x.dist <- x.dist
  res$x.class2 <- x.class2
  res$hull.all <- hull.all
  res$FRic.all <- FRic.all
  if(exists("warning", where = .GlobalEnv) && !is.function(warning)){
    res$warning <- warning
  } else {
    res$warning <- FALSE
  }
  
  
  return(res)
}