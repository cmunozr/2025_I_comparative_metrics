#' @title Pre-compute Stand-Level Intrinsic Metrics
#' @description Computes single-stand richness, Rao's Quadratic Entropy, and Functional Richness
#' across all posterior draws prior to pairwise comparisons.
precalculate_stand_metrics <- function(pred_mat, traits, divmax) {
  # pred_mat: [n_posteriors x n_species]
  num_posteriors <- nrow(pred_mat)
  
  mask <- pred_mat > 0
  prob <- 1 - exp(-pred_mat)
  
  # Unforced expected richness per posterior
  es_unforced <- rowSums(prob * mask)
  
  # Rao's Quadratic Entropy per posterior
  rao_q <- divc_calc(mask, tr = traits, matrix.ver = TRUE, scalar = divmax)
  
  # Functional Richness per posterior
  fric <- numeric(num_posteriors)
  for (p in seq_len(num_posteriors)) {
    fric[p] <- fric_process_row(mask[p, ], traits)
  }
  
  list(
    pred_mat = pred_mat,
    mask = mask,
    prob = prob,
    es_unforced = es_unforced,
    rao_q = rao_q,
    fric = fric
  )
}

#' @title Pairwise Evaluation
#' @description Evaluates interaction metrics and combines them with pre-calculated metrics.
calculate_pairwise_metrics_fast <- function(metso_prep, bau_prep, num_posteriors, xi = 1e-3) {
  
  # 1. Scalar differences from pre-computed metrics
  E_PDF_unforced <- 1 - (bau_prep$es_unforced / metso_prep$es_unforced)
  E_RaoQ_loss    <- 1 - (bau_prep$rao_q / metso_prep$rao_q) 
  
  # Handle FRic loss safely
  fric_m <- metso_prep$fric
  fric_b <- bau_prep$fric
  valid_fric <- !is.na(fric_m) & fric_m > 0
  E_FRic_loss <- rep(NA_real_, num_posteriors)
  E_FRic_loss[valid_fric] <- 1 - (fric_b[valid_fric] / fric_m[valid_fric])
  
  # 2. Pairwise interaction metrics
  E_BCD <- bray_curtis(metso_prep$pred_mat, bau_prep$pred_mat)
  
  E_PDF_forced <- numeric(num_posteriors)
  E_MSA_loss   <- numeric(num_posteriors)
  EG_loss      <- numeric(num_posteriors)
  
  for (p in seq_len(num_posteriors)) {
    m_metso <- metso_prep$mask[p, ]
    es_m    <- metso_prep$es_unforced[p]
    
    if (is.na(es_m) || es_m == 0) {
      E_PDF_forced[p] <- NA
      E_MSA_loss[p]   <- NA
      EG_loss[p]      <- NA
      next
    }
    
    # Forced expected richness of BAU under METSO's present species
    es_bau_forced <- sum(bau_prep$prob[p, m_metso])
    E_PDF_forced[p] <- 1 - (es_bau_forced / es_m)
    
    # MSA Loss
    e_m <- metso_prep$pred_mat[p, m_metso]
    e_b <- bau_prep$pred_mat[p, m_metso]
    ep_m <- metso_prep$prob[p, m_metso]
    
    ratio_abun <- e_b / e_m
    E_MSA_loss[p] <- 1 - (sum(ep_m * pmin(ratio_abun, 1)) / es_m)
    
    # EG Loss
    log_ratio <- log((e_b + xi) / (e_m + xi))
    EG_loss[p] <- 1 - exp(sum(log_ratio) / es_m)
  }
  
  data.frame(
    posterior = rep(seq_len(num_posteriors), each = 7),
    metrics = rep(c("E_PDF_forced", "E_PDF_unforced", "E_MSA_loss", "EG_loss", 
                    "E_FRic_loss", "E_RaoQ_loss", "E_BCD"), times = num_posteriors),
    val = c(rbind(E_PDF_forced, E_PDF_unforced, E_MSA_loss, EG_loss, 
                  E_FRic_loss, E_RaoQ_loss, E_BCD))
  )
}

#' @title Vectorized Calculation of Expected Biodiversity Metrics
#' @description Computes an expected suite of taxonomic and functional biodiversity metrics comparing business-as-usual (BAU) versus METSO baseline scenarios.
#' Metric formulations are mapped to the methodological equations in the manuscript.
#' @param predY_metso_mat Matrix of posterior predictions for the natural/baseline scenario.
#' @param predY_bau_mat Matrix of posterior predictions for the pressure/management scenario.
#' @param Traits Processed traits dataset from dbFD_preprocess_traits.
#' @param xi Stabilizing constant used for Expected Geometric Mean Abundance calculations to prevent mathematical instability at zero.
#' @return A data frame in long format recording metrics per posterior draw.

calculate_metrics_vectorized <- function(predY_metso_mat, predY_bau_mat, Traits = TrData_processed, xi = 1e-3) {
  
  if(!is.matrix(predY_metso_mat)) predY_metso_mat <- as.matrix(predY_metso_mat)
  if(!is.matrix(predY_bau_mat)) predY_bau_mat <- as.matrix(predY_bau_mat) 
  
  num_posteriors <- nrow(predY_metso_mat)
  
  # Pre-allocate results
  E_PDF_forced <- numeric(num_posteriors)
  E_PDF_unforced <- numeric(num_posteriors)
  E_MSA_loss <- numeric(num_posteriors)
  EG_loss <- numeric(num_posteriors)
  E_FRic_loss <- numeric(num_posteriors)
  
  # 1 & 2: Define masks (alpha is incorporated already in the predictions)
  mask_metso <- predY_metso_mat > 0
  mask_bau <- predY_bau_mat > 0
  
  # 3: Calculate expected probabilities for richness sums
  prob_metso <- 1 - exp(-predY_metso_mat)
  prob_bau <- 1 - exp(-predY_bau_mat)
  
  # Some metrics can be calculated at the complete posterior level, 
  # rather than per chunk of posterior, for efficiency.

  # Expected Loss of Rao's Quadratic Entropy (Eq. 22 and 23)
  divmax <- ade4::divcmax(as.dist(Traits$x.dist))$value
  E_RaoQ_loss <- 1 - (divc_calc(mask_bau, tr = Traits, scalar = divmax) / divc_calc(mask_metso, tr = Traits, scalar = divmax))
  
  # Expected Bray-Curtis Dissimilarity (1 - Sørensen similarity) (Eq. 12, 13, 14)
  E_BCD <- bray_curtis(predY_metso_mat, predY_bau_mat)
  
  # Process each posterior
  for(p in seq_len(num_posteriors)) {
    m_metso <- mask_metso[p, ]
    m_bau <- mask_bau[p, ]
    
    # Expected Richness Extraction
    ES_metso <- sum(prob_metso[p, m_metso])
    ES_bau_unforced <- sum(prob_bau[p, m_bau])
    ES_bau_forced <- sum(prob_bau[p, m_metso])
    
    # Mathematical safeguard: Skip iteration if baseline richness is zero
    if (ES_metso == 0) {
      E_PDF_forced[p] <- NA
      E_PDF_unforced[p] <- NA
      E_MSA_loss[p] <- NA
      EG_loss[p] <- NA
      E_FRic_loss[p] <- NA
      next
    }
    
    # Potentially Disappeared Fraction (Eq. 4 and 5)
    E_PDF_forced[p] <- 1 - (ES_bau_forced / ES_metso)
    E_PDF_unforced[p] <- 1 - (ES_bau_unforced / ES_metso)
    
    # Extract expected abundances for species present in baseline
    E_metso <- predY_metso_mat[p, m_metso]
    E_bau_forced_abun <- predY_bau_mat[p, m_metso]
    EP_metso <- prob_metso[p, m_metso]
    
    # Mean Species Abundance Loss (Eq. 8)
    ratio_abun <- E_bau_forced_abun / E_metso
    E_MSA_loss[p] <- 1 - (sum(EP_metso * pmin(ratio_abun, 1)) / ES_metso)
    
    # Expected Geometric Mean Abundance Loss (Eq. 10 and 11)
    # Restricted to the baseline community and utilizes stabilizing constant (xi)
    log_ratio_abun <- log((E_bau_forced_abun + xi) / (E_metso + xi))
    EG_j <- exp(sum(log_ratio_abun) / ES_metso)
    EG_loss[p] <- 1 - EG_j
    
    # Expected Functional Richness Loss (Eq. 16 and 17)
    FRic_metso <- fric_process_row(m_metso, Traits)
    FRic_bau <- fric_process_row(m_bau, Traits)
    
    if (is.na(FRic_metso) || FRic_metso == 0) {
      E_FRic_loss[p] <- NA
    } else {
      E_FRic_loss[p] <- 1 - (FRic_bau / FRic_metso)
    }
  }
  
  data.frame(
    posterior = rep(seq_len(num_posteriors), each = 7),
    metrics = rep(c("E_PDF_forced", "E_PDF_unforced", "E_MSA_loss", "EG_loss", "E_FRic_loss", "E_RaoQ_loss", "E_BCD"), 
                  times = num_posteriors),
    val = c(rbind(E_PDF_forced, E_PDF_unforced, E_MSA_loss, EG_loss, E_FRic_loss, E_RaoQ_loss, E_BCD))
  )
}

#-----------------

#' @title Expected Bray-Curtis dissimilarity (1 - quantitative Sørensen similarity)
#' @description Computes the Bray-Curtis dissimilarity between baseline and scenario expected abundances.
#' Mathematically equivalent to 1 - quantitative Sørensen similarity.
#' See Equations 12, 13, and 14 in the manuscript.
#' @param baseline_mat Numeric matrix of expected abundances for the baseline scenario.
#' @param scenario_mat Numeric matrix of expected abundances for the pressure scenario.
#' @param array Logical; whether the inputs are 3D arrays. Defaults to FALSE.
#' @return A numeric vector of Bray-Curtis dissimilarity values.
bray_curtis <- function(baseline_mat, scenario_mat, array = FALSE) {
  
  if (!identical(dim(baseline_mat), dim(scenario_mat))) {
    stop("Dimensions of 'baseline_mat' and 'scenario_mat' must be identical.")
  }
  
  # Eq. 12 and 14: sum(|scenario - baseline|) / sum(scenario + baseline)
  # Mathematically equivalent to 1 - quantitative Sørensen (Eq. 13)
  sum_abs_diff <- abs(scenario_mat - baseline_mat)
  total_abundance <- abs(scenario_mat + baseline_mat)
  
  if (array) {
    total_abundance <- apply(total_abundance, 3, rowSums)
    sum_abs_diff <- apply(sum_abs_diff, 3, rowSums)
  } else {
    total_abundance <- rowSums(total_abundance, dims = 1)
    sum_abs_diff <- rowSums(sum_abs_diff, dims = 1)
  }
  
  # Calculate dissimilarity ratio
  bray_curtis_dissimilarity <- sum_abs_diff / total_abundance
  
  # Handle Inf and NaN for empty or zero-abundance communities
  bray_curtis_dissimilarity <- ifelse(
    is.infinite(bray_curtis_dissimilarity) | is.nan(bray_curtis_dissimilarity), 
    NA, 
    bray_curtis_dissimilarity
  )
  
  return(bray_curtis_dissimilarity)
}

#------------------

#' @title Functional Richness Processing
#' @description Calculates the functional richness (convex hull volume) for a single site/row.
#' Based on https://doi.org/10.1890/08-2244.1. Taken from ade4.
#' @param site.row Logical or numeric vector indicating species presence at a site.
#' @param tr Processed traits object containing 'traits.FRic', 'FRic.all', etc.
#' @param std.fric Logical; whether to standardize the functional richness. Defaults to TRUE.
#' @return Numeric value representing the functional richness volume.
fric_process_row <- function(site.row, tr, std.fric = T) {
  
  traits <- tr$traits.FRic
  x.class2 <- tr$x.class2
  hull.all <- tr$hull.all
  FRic.all <- tr$FRic.all
  war <- tr$warning
  
  sp_indices <- which(site.row > 0)
  nb.sp <- length(sp_indices)
  
  # Subset traits for this site. Note: drop=FALSE ensures it stays a matrix even if only 1 species is present
  tr.FRic <- traits[sp_indices, , drop = FALSE]
  
  # Initialize return value
  FRic_value <- NA
  
  # Check Minimum Species Count
  if (nb.sp < 3) {
    return(NA)
  }
  
  # Factor / Ordered Traits
  if (all(x.class2 == "factor" | x.class2 == "ordered")) {
    
    if (length(x.class2) == 1 & x.class2[1] == "ordered") {
      # Single Ordered Trait
      tr.range <- range(tr.FRic[, 1])
      t.range <- tr.range[2] - tr.range[1]
      
      if (!std.fric) FRic_value <- t.range
      if (std.fric)  FRic_value <- t.range / FRic.all
      
    } else {
      # Categorical / Mixed
      if (!std.fric) FRic_value <- nrow(unique(tr.FRic))
      if (std.fric)  FRic_value <- nrow(unique(tr.FRic)) / FRic.all
    }
    
  } else {
    
    # Numeric / Continuous Traits
    
    # CASE: Multi-dimensional (>1 trait)
    if (dim(tr.FRic)[2] > 1 & nb.sp >= 3) {
      
      # Threshold logic from your original code
      if (war)  thresh <- 4
      if (!war) thresh <- 3
      
      if (nb.sp >= thresh) {
        
        # tryCatch Block: Catches geometry errors (e.g., if species are collinear/coplanar)
        hull_vol <- tryCatch({
          geometry::convhulln(tr.FRic, "FA")$vol
        }, error = function(e) {
          return(NA)
        })
        
        # assign if calculation succeeded
        if (!is.na(hull_vol)) {
          if (!std.fric) FRic_value <- hull_vol
          if (std.fric)  FRic_value <- hull_vol / FRic.all
        }
      }
    }
    
    # Single dimension (1 trait)
    if (dim(tr.FRic)[2] == 1) {
      tr.range <- range(tr.FRic[, 1])
      t.range <- tr.range[2] - tr.range[1]
      
      if (!std.fric) FRic_value <- t.range
      if (std.fric)  FRic_value <- t.range / FRic.all
    }
  }
  
  return(FRic_value)
}

#-----------------------------

#' @title Rao's Quadratic Entropy Coefficient Calculation
#' @description Represents the average distance between two randomly selected individuals from the community.
#' Based on https://doi.org/10.1890/08-2244.1. Taken from ade4.
#' @param mat Matrix of expected abundances or presence masks.
#' @param tr Processed traits object containing distance matrix 'x.dist'.
#' @param matrix.ver Logical; use fast matrix algebra path. Defaults to TRUE.
#' @param scalar Numeric value to standardize RaoQ against the theoretical maximum.
#' @return Numeric vector of RaoQ values per site.
divc_calc <- function(mat, tr, matrix.ver = TRUE, scalar) {
  
  # Setup Distance Matrix
  # Check if distance exists in the trait object
  if (is.null(tr$x.dist)) {
    # Fallback: If no distance provided, assume equidistance (Gini-Simpson equivalent)
    # Create a distance matrix where diag=0 and off-diag=sqrt(2)
    n_sp <- ncol(mat)
    dis_matrix <- (matrix(1, n_sp, n_sp) - diag(rep(1, n_sp))) * sqrt(2)
  } else {
    if (!inherits(tr$x.dist, "dist")) stop("Object of class 'dist' expected for distance")
    dis_matrix <- as.matrix(tr$x.dist)
  }
  
  # Input Integrity Checks
  if (any(mat < 0, na.rm = TRUE)) stop("Negative value in abundance matrix")
  
  # matrix.ver: Matrix Algebra (The "Fast" Path)
  if (matrix.ver) {
    
    D_sq <- dis_matrix^2
    
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    
    S <- rowSums(mat)
    
    # Handle sites with 0 abundance 
    S[S == 0] <- NA
    
    # The Vectorized Calculation:
    # No looping per every site, calculate the diagonal 
    # A. (mat %*% dis) multiplies richness by Distance
    # B. * mat multiplies that result by richness element-wise
    # C. rowSums aggregates it to get the Rao value per site
    
    numerator <- rowSums((mat %*% D_sq) * mat)
    
    results <- numerator / (2 * (S^2))
    
  } else {
    
    # Data Frame (The "Safe/Legacy" Path)
    # Transpose to Species x Sites 
    df <- as.data.frame(t(mat))
    dis_dist <- as.dist(dis_matrix)
    
    results <- sapply(1:ncol(df), function(i) {
      if (sum(df[, i]) < 1e-16) {
        return(0)
      } else {
        # Original single-site formula
        val <- (t(df[, i]) %*% (as.matrix(dis_dist)^2) %*% df[, i]) / (2 * (sum(df[, i])^2))
        return(as.numeric(val))
      }
    })
  }
  
  if(!is.null(scalar)){
    results <- results/scalar
  }
  
  return(results)
}

#-----------------

#' @title Community Weighted Mean (Parallelized)
#' @description Calculates trait means weighted by abundance. Note: Considered likely unfeasible for current framework evaluation.
#' Based on https://doi.org/10.1890/08-2244.1
#' @param pred.object Array of predictions.
#' @param trait.processed.object Processed traits data.
#' @param cwm.type Character string for CWM type ("dom" or "all").
#' @param parallel Logical; whether to use multicore processing.
#' @param use.cores Integer; number of cores to use.
#' @return List of weighted traits.
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

# utilities

#' @title Safe Logarithm for Zero Values
#' @description Applies a log transformation while safely maintaining exact zeros.
#' @param x Numeric vector.
#' @return Log-transformed vector with zeros preserved.
log_zero <- function(x) {
  ifelse(x == 0, 0, log(x))
}

#----------------------------

# preprocess traits. Taken from FD function dbFD (original function quite slowly)
##  https://doi.org/10.1890/08-2244.1

#' @title Preprocess Functional Traits (Optimized dbFD)
#' @description A helper function derived from FD::dbFD that pre-processes functional trait data.
#' It computes the species-by-species distance matrix, performs PCoA to build the functional space, 
#' and calculates the global convex hull volume. Designed for efficiency in simulation or Bayesian workflows.
#' @param x Trait data.
#' @param a Abundance matrix (dummy matrix used for setup).
#' @param w Trait weights.
#' @param w.abun Logical.
#' @param stand.x Logical; standardize traits.
#' @param ord Ordination method.
#' @param asym.bin Asymmetric binary handling.
#' @param corr Distance correction method.
#' @param calc.FRic Logical; calculate functional richness bounds.
#' @param m Dimensions for PCoA.
#' @param stand.FRic Logical.
#' @param scale.RaoQ Logical.
#' @param calc.FGR Logical.
#' @param clust.type Clustering method.
#' @param km.inf.gr K-means lower bound.
#' @param km.sup.gr K-means upper bound.
#' @param km.iter K-means iterations.
#' @param km.crit K-means criteria.
#' @param calc.CWM Logical.
#' @param CWM.type Character.
#' @param calc.FDiv Logical.
#' @param dist.bin Distance binary parameter.
#' @param print.pco Logical.
#' @param messages Logical; print console warnings.
#' @return A list of preprocessed 'ingredients' for downstream diversity metric calculation.
dbFD_preprocess_traits <- function (x, a, w, w.abun = TRUE, stand.x = TRUE, ord = c("podani", "metric"), 
                                    asym.bin = NULL, corr = "sqrt", calc.FRic = TRUE, m = "max", stand.FRic = FALSE, 
                                    scale.RaoQ = FALSE, calc.FGR = F, clust.type = "ward", 
                                    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, 
                                    km.crit = c("calinski", "ssi"), calc.CWM = TRUE, CWM.type = c("dom", "all"), 
                                    calc.FDiv = F, dist.bin = 2, print.pco = FALSE, messages = TRUE) 
{
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
      stop("Different number of species in 'x' and 'a'.", "\n")
    if (any(ab.t.row != x.rn)) 
      stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", "\n")
  }
  
  a <- as.matrix(a)
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    stop("At least one community has zero-sum abundances (no species).", "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    stop("At least one species does not occur in any community (zero total abundance across all communities).", "\n")
  if (!missing(w) & is.dist.x) 
    stop("When 'x' is a distance matrix, 'w' should be left missing.", "\n")
  if (!missing(w) & !is.dist.x) {
    if (!is.numeric(w) | length(w) != t.x) 
      stop("'w' should be a numeric vector of length = number of traits.", "\n")
    else w <- w/sum(w)
  }
  if (missing(w)) 
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) 
        x[, x.class == "character"] <- as.factor(x[, x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        } else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      } else {
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
            cat("Warning: Species with missing trait values have been excluded.", "\n")
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
            cat("Warning: Species with missing trait values have been excluded.", "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        } else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)])) 
            stop("'dist.bin' must be an integer between 1 and 10.", "\n")
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
        cat("Warning: Species with missing trait values have been excluded.", "\n")
    } else x <- x
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
        cat("Warning: Species with missing trait values have been excluded.", "\n")
    } else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.", "\n")
    } else x <- x
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
        cat("Warning: Species with missing trait values have been excluded.", "\n")
    } else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x))) 
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).", "\n")
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
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2)) 
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", "\n")
      if (is.euclid(x.dist2)) 
        if (messages) 
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages) 
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", "\n")
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
      cat("FEVe: Could not be calculated for communities with <3 functionally singular species.", "\n")
  if (min.nb.sp < 2) 
    if (messages) 
      cat("FDis: Equals 0 in communities with only one functionally singular species.", "\n")
  if (calc.FRic) {
    x.class2 <- sapply(x, data.class)
    if (all(x.class2 == "factor" | x.class2 == "ordered")) {
      if (length(x.class2) == 1 & x.class2[1] == "ordered") {
        traits.FRic1 <- rank(x[, 1])
        names(traits.FRic1) <- x.rn
        traits.FRic <- data.frame(traits.FRic1)
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one ordinal trait present in 'x'. FRic was measured as the range of the ranks, NOT as the convex hull volume.", "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when 'x' is a single ordinal trait.", "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      } else {
        traits.FRic <- x
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only categorical and/or ordinal trait(s) present in 'x'. FRic was measured as the number of unique trait combinations, NOT as the convex hull volume.", "\n")
        if (stand.FRic) 
          FRic.all <- nrow((unique(traits.FRic)))
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when only categorical and/or ordinal trait(s) present in 'x'.", "\n")
        }
      }
    } else {
      if (x.pco$nf == 1) {
        traits.FRic <- x.pco$li
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.", "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.", "\n")
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
              cat("FRic: To respect s >= 2^t, FRic could not be calculated for communities with <4 functionally singular species.", "\n")
          }
          else m.min <- floor(log2(min.nb.sp))
        } else {
          if (min.nb.sp < 3) {
            nb.sp2 <- nb.sp[nb.sp > 2]
            m.max <- min(nb.sp2) - 1
            if (messages) 
              cat("FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species.", "\n")
          }
          else m.max <- m.max
        }
        if (is.numeric(m) & m <= 1) 
          stop("When 'm' is an integer, it must be >1.", "\n")
        if (is.numeric(m) & m > m.max) 
          m <- m.max
        if (m == "min") 
          m <- m.min
        if (m == "max") 
          m <- m.max
        if (!is.numeric(m) & m != "min" & m != "max") 
          stop("'m' must be an integer >1, 'min', or 'max'.", "\n")
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
              cat("FRic: Quality of the reduced-space representation =", qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr != "none") {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (based on corrected distance matrix) =", qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr == "none") {
            delta <- -0.5 * bicenter.wt(x.dist * x.dist)
            lambda <- eigen(delta, symmetric = TRUE, only.values = TRUE)$values
            sum.m <- sum(lambda[1:m])
            sum.n <- sum(lambda)
            lambda.neg <- c(lambda[lambda < 0])
            max.neg <- abs(min(lambda.neg))
            qual.FRic <- (sum.m + (length(lambda[1:m]) * max.neg))/(sum.n + ((length(lambda) - 1) * max.neg))
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) =", qual.FRic, "\n")
          }
        }
        if (m >= x.pco$nf) {
          qual.FRic = 1
          traits.FRic <- x.pco$li
          if (x.pco$nf == 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. The 2 PCoA axes were kept as 'traits'.", "\n")
          if (x.pco$nf > 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. All", x.pco$nf, "PCoA axes were kept as 'traits'.", "\n")
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
  res$eigenvalues <- x.pco$eig
  if(exists("warning", where = .GlobalEnv) && !is.function(warning)){
    res$warning <- warning
  } else {
    res$warning <- FALSE
  }
  
  return(res)
}

#' Get or Create Test Set by Sampling Complete Strata Groups
#'
#' @param gpkg_path Character. Path to the full stand GPKG layer.
#' @param test_group_path Character. Path to RDS where selected stand metadata is cached.
#' @param n_groups Integer. Number of groups to sample (default = 20).
#' @param min_stands Integer. Minimum stands per group per scenario (default = 10).
#' @param max_stands Integer. Maximum stands per group per scenario (default = 20).
#' @param seed Integer. Random seed for reproducible group sampling.
#'
#' @return A data frame containing standid, metso, and group_number for the sampled groups.
get_or_create_test_groups <- function(gpkg_path,
                                      test_group_path = file.path("results", "test_group_stands.rds"),
                                      n_groups = 20,
                                      min_stands = 10,
                                      max_stands = 20,
                                      seed = 11072024) {
  
  if (file.exists(test_group_path)) {
    message("Loading existing test group set from: ", test_group_path)
    return(readRDS(test_group_path))
  }
  
  message("Generating new balanced test group set with seed ", seed, "...")
  set.seed(seed)
  
  # Read full layer (both METSO and control stands)
  sp_all <- sf::read_sf(gpkg_path) |>
    sf::st_drop_geometry()
  
  # Build identical group indices used in S08b
  grouped_stands <- sp_all |>
    dplyr::group_by(year, regional_group, treespecies) |>
    dplyr::mutate(group_number = dplyr::cur_group_id()) |>
    dplyr::ungroup()
  
  # Count stands available per group separately for METSO and control
  group_counts <- grouped_stands |>
    dplyr::group_by(group_number, metso) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = metso, values_from = n, values_fill = 0, names_prefix = "metso_")
  
  # Filter groups meeting criteria in both scenarios
  # metso_1 = METSO, metso_0 = Control
  eligible_groups <- group_counts |>
    dplyr::filter(
      metso_1 >= min_stands & metso_1 <= max_stands,
      metso_0 >= min_stands & metso_0 <= max_stands
    ) |>
    dplyr::pull(group_number)
  
  if (length(eligible_groups) == 0) {
    stop("No groups found satisfying stand thresholds between ", min_stands, " and ", max_stands, 
         " for both treatment and control.")
  }
  
  n_sample <- min(n_groups, length(eligible_groups))
  selected_groups <- sample(eligible_groups, size = n_sample)
  
  # Isolate all stands belonging to the selected groups
  test_selection <- grouped_stands |>
    dplyr::filter(group_number %in% selected_groups) |>
    dplyr::select(standid, metso, group_number, year, regional_group, treespecies)
  
  # Ensure target directory exists
  out_dir <- dirname(test_group_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  saveRDS(test_selection, file = test_group_path)
  message(sprintf("Saved %d stands across %d groups to %s", nrow(test_selection), n_sample, test_group_path))
  
  return(test_selection)
}