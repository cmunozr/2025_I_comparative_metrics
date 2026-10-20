#' DESCRIPTION:
#' S10a: Evaluation of METSO conservation impact across biodiversity metrics
#' using fixed-effect Ordinary Least Squares (fixest::feols).
#' 
#' METHODOLOGY:
#' - Evaluates whether the conditional mean difference (val_mean) between METSO 
#'   and matched BAU controls significantly differs from zero across strata.
#' - Specification: val_mean ~ 1 + i(reg) + i(trees) + i(year), cluster = ~group_number
#' - Reference Test: Evaluates intercept on posterior mean metrics.
#' - Robustness: Propagates parameter and prediction uncertainty across posterior draws.
#' - Draw-level national and region-standardized marginal means from the same weighted OLS.
#' 
#' INPUTS:
#' - Group Parquet files: results/metrics/[expected_type]/[run_mode]/metrics_group_*.parquet
#' - Group Matches: 
#'  - results/matched_pairs_groups_metso.rds
#'  - results/matched_pairs_groups_control.rds
#' 
#' OUTPUTS:
#' - Summary Table: results/causal_analysis/rq1_causal_effects_summary.rds
#' - Posterior Draws: results/causal_analysis/rq1_draws_robustness.rds
#' - Regional Draws / Summary: results/causal_analysis/rq1_draws_regional_*.rds

library(tidyverse)
library(arrow)
library(data.table)
library(fixest)
library(here)

source(file.path("code", "config_model.R"))

cat(Sys.time())

# --- 1. Environment & Directory Setup ---
test <- FALSE
expected_string <- "_expected_true"

os <- Sys.info()["sysname"]
if (os == "Windows") {
  wr_dir <- "results"
} else if (os %in% c("Linux", "Darwin")) {
  node <- Sys.info()["nodename"]
  if (node == "lema") {
    wr_dir <- "results"
  } else if (node == "gurnah") {
    wr_dir <- "/export/scratch/tmp/munozcs/results"
  } else {
    wr_dir <- "results"
  }
}

run_mode <- "draws" # flag to remind

if (run_mode == "mean") {
  stop("Legacy 'mean' mode is disabled for production to avoid metric bias (e.g. FRic). Use 'draws' or 'all'.")
  
  mean_dir <- if (test) {
    file.path(wr_dir, "metrics", "test", expected_string, "mean")
  } else {
    file.path(wr_dir, "metrics", expected_string, "mean")
  }
}


draws_dir <- if (test) {
  file.path(wr_dir, "metrics", "test", expected_string, "draws")
} else {
  file.path(wr_dir, "metrics", expected_string, "draws")
}

out_causal_dir <- file.path(wr_dir, "causal_analysis")
dir.create(out_causal_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load Stratum Metadata for Fixed Factors ---
# S09 parquet outputs contain standid_treated, group_number, posterior, metrics, val_mean, val_sd, n_controls.
# Join the matched group covariates (reg, trees, year) for fixed effect controls.

matches_metso <- readRDS(file.path(wr_dir, "matched_pairs_groups_metso.rds"))
matches_control <- readRDS(file.path(wr_dir, "matched_pairs_groups_control.rds"))
groups_meta <- unique(
  bind_rows(matches_metso, matches_control)[, c("group_number", "reg", "trees", "year")]
)
setDT(groups_meta)
groups_meta[, `:=`(
  reg   = as.factor(reg),
  trees = as.factor(trees),
  year  = as.factor(year)
)]

# --- 3. Mode Mean ---
cat("\n--- Evaluating Posterior Mean Metrics ---\n")

if(run_mode == "mean"){
  
  ds_mean <- arrow::open_dataset(mean_dir)
  
  # Collect data
  dt_mean <- ds_mean |>
    select(standid_treated, group_number, posterior, metrics, val_mean, val_sd, n_controls) |>
    mutate(weigth_sd = (1/val_sd)^2) |> 
    collect()
  setDT(dt_mean)
  
  dt_mean <- merge(dt_mean, groups_meta, by = "group_number", all.x = TRUE)
  
  metrics_list <- sort(unique(dt_mean$metrics))
  
  results_q1 <- list()
  q1_mean_fitted_models <- list()
  
  for (m in metrics_list) {
    # m <- metrics_list[1]  # For testing/debugging
    sub_df <- dt_mean[metrics == m]
    
    if (nrow(sub_df) == 0) next
    
    # Linear model with environmental strata fixed factors and clustered SEs
    fit <- tryCatch({
      feols(
        val_mean ~ 1 + i(reg) + i(trees) + i(year),
        data    = sub_df,
        cluster = ~group_number#, 
        #weights = ~weigth_sd 
      )
    }, error = function(e) {
      message("Model failed for metric ", m, ": ", e$message)
      return(NULL)
    })
    
    q1_mean_fitted_models[[m]] <- fit
    
    if (!is.null(fit)) {
      
      int_coef <- as.numeric(coef(fit)["(Intercept)"])
      int_se   <- as.numeric(se(fit)["(Intercept)"])
      int_pval <- as.numeric(pvalue(fit)["(Intercept)"])
      ci_vals  <- as.numeric(confint(fit)["(Intercept)", ])
      r2_val   <- as.numeric(r2(fit, "ar2"))
      
      row_dt <- data.table(
        metric       = as.character(m),
        n_stands     = as.integer(nrow(sub_df)),
        n_groups     = as.integer(uniqueN(sub_df$group_number)),
        intercept    = int_coef,
        se           = int_se,
        t_stat       = int_coef / int_se,
        p_value      = int_pval,
        ci_2.5       = ci_vals[1],
        ci_97.5      = ci_vals[2],
        r2_adjusted  = r2_val
      )
      
      results_q1[[m]] <- row_dt
    }
  }
  
  q1_mean_summary <- rbindlist(results_q1)
  
  saveRDS(q1_mean_summary, file.path(out_causal_dir, "q1_causal_effects_summary.rds"))
  saveRDS(q1_mean_fitted_models, file.path(out_causal_dir, "q1_fitted_models.rds"))
}

# --- 4. Mode Draws ---
cat("\n--- Evaluating Posterior Draws Uncertainty ---\n")

ds_draws <- arrow::open_dataset(draws_dir)

dt_draws <- ds_draws |>
  select(standid_treated, group_number, posterior, metrics, val_mean, val_sd) |>
  mutate(weigth_sd = (1/val_sd)^2) |>
  collect()
setDT(dt_draws)

metrics_list <- sort(unique(dt_draws$metrics))

dt_draws <- merge(dt_draws, groups_meta, by = "group_number", all.x = TRUE)

draw_ids <- sort(unique(dt_draws$posterior))

# Target population: unique matched METSO stands. Marginal averaging uses equal
# stand weights, independently of the inverse-variance weights used to fit OLS.
stand_pop <- unique(dt_draws[, .(standid_treated, group_number, reg, trees, year)])

national_grid <- stand_pop[, .(n_target = .N), by = .(reg, trees, year)]
national_grid[, share := n_target / sum(n_target)]

tree_year_mix <- stand_pop[, .(n_target = .N), by = .(trees, year)]
tree_year_mix[, share := n_target / sum(n_target)]

regional_grid <- as.data.table(merge(
  data.frame(reg = unique(stand_pop$reg)),
  as.data.frame(tree_year_mix), by = NULL
))
regional_counts <- stand_pop[, .(
  n_region_stands = .N,
  n_region_groups = uniqueN(group_number)
), by = reg]

# A marginal mean is a linear combination of the fitted OLS coefficients.
# Its SE uses the full matched-group clustered covariance, not mean(prediction SEs).

marginal_estimate <- function(fit, grid, weight_col = "share") {
  
  w <- grid[[weight_col]]
  w <- w / sum(w)
  
  pred <- predict(fit, newdata = as.data.frame(grid))
  estimate <- sum(w * pred)
  
  b <- coef(fit)
  V <- vcov(fit)
  
  x <- matrix(
    0,
    nrow = nrow(grid),
    ncol = length(b),
    dimnames = list(NULL, names(b))
  )
  
  x[, "(Intercept)"] <- 1
  
  for (nm in names(b)) {
    
    if (startsWith(nm, "reg::")) {
      level <- sub("reg::", "", nm, fixed = TRUE)
      x[, nm] <- as.character(grid$reg) == level
    }
    
    if (startsWith(nm, "trees::")) {
      level <- sub("trees::", "", nm, fixed = TRUE)
      x[, nm] <- as.character(grid$trees) == level
    }
    
    if (startsWith(nm, "year::")) {
      level <- sub("year::", "", nm, fixed = TRUE)
      x[, nm] <- as.character(grid$year) == level
    }
  }
  
  # Average design vector for the requested target population.
  l <- colSums(x * w)
  
  # Safety check: our coefficient mapping must reproduce fixest predictions.
  estimate_from_coef <- sum(l * b)
  
  if (!isTRUE(all.equal(
    estimate,
    estimate_from_coef,
    tolerance = 1e-10
  ))) {
    stop("Marginal coefficient mapping does not reproduce fixest predictions.")
  }
  
  se <- sqrt(drop(t(l) %*% V %*% l))
  
  tibble(
    estimate = estimate,
    se = se,
    ci_2.5 = estimate - 1.96 * se,
    ci_97.5 = estimate + 1.96 * se
  )
}

# Pre-allocate
total_iterations <- length(metrics_list) * length(draw_ids)
draws_results <- vector("list", total_iterations)
regional_results <- vector("list", total_iterations)
counter <- 1

for (m in metrics_list) {
  # m <- "E_BCD"
  sub_m <- dt_draws[metrics == m]
  
  for (d in draw_ids) {
    # d <- 11
    sub_md <- sub_m[posterior == d]
    if (nrow(sub_md) == 0) next
    
    fit_d <- tryCatch({
      feols(
        val_mean ~ 1 + i(reg) + i(trees) + i(year),
        data    = sub_md,
        cluster = ~group_number, 
        weights = ~weigth_sd
      )
    }, error = function(e) {
      message(sprintf("Model failed for metric %s (draw %s): %s", m, d, e$message))
      return(NULL)
    })
    
    if (!is.null(fit_d)) {
      # National: average fitted contrasts over the observed METSO stand mix.
      national_margin <- marginal_estimate(fit_d, national_grid)
      row_draw <- cbind(data.table(
        metric    = as.character(m),
        posterior = as.integer(d),
        n_stands  = nrow(stand_pop),
        n_groups  = uniqueN(stand_pop$group_number),
        n_fitted  = nobs(fit_d)
      ), national_margin)
      
      # Regional: set region, then average over the SAME national joint
      # tree-type/year mix in every region (not each region's own composition).
      # Unobserved region x tree/year cells rely on the additive OLS assumption.
      region_margin <- rbindlist(lapply(split(regional_grid, regional_grid$reg, drop = TRUE), function(g) {
        cbind(data.table(reg = g$reg[1]), marginal_estimate(fit = fit_d, grid = g))
      }))
      region_margin <- merge(region_margin, regional_counts, by = "reg")
      
      region_margin <- region_margin |> 
        rename(marginal_mean = estimate) |>
        mutate(
          metric = as.character(m),
          posterior = as.integer(d),
          t_stat = marginal_mean / se,
          p_value = 2 * pnorm(abs(t_stat), lower.tail = FALSE)
        )
      setcolorder(region_margin, c(
        "metric", "posterior", "reg", "n_region_stands", "n_region_groups",
        "marginal_mean", "se", "t_stat", "p_value", "ci_2.5", "ci_97.5"
      ))
      
      draws_results[[counter]] <- row_draw
      regional_results[[counter]] <- region_margin
      
      counter <- counter + 1
    }
  }
}

# Remove unused pre-allocated slots and bind
draws_dt_all <- rbindlist(draws_results[seq_len(counter - 1)])
regional_dt_all <- rbindlist(regional_results[seq_len(counter - 1)])

# Across-draw intervals below describe JSDM posterior uncertainty in marginal
# point estimates. Per-draw ci_2.5/ci_97.5 separately use clustered OLS SEs.
q1_robustness <- draws_dt_all[, .(
  n_draws_evaluated    = .N,
  estimate_mean   = mean(estimate, na.rm = TRUE),
  estimate_median = median(estimate, na.rm = TRUE),
  sd_estimate     = sd(estimate, na.rm = TRUE),
  ci_2.5               = quantile(estimate, 0.025, na.rm = TRUE),
  ci_97.5              = quantile(estimate, 0.975, na.rm = TRUE),
  prob_positive        = mean(estimate > 0, na.rm = TRUE),
  sign_stable          = (all(estimate > 0) || all(estimate < 0))
), by = .(metric)]

regional_summary <- regional_dt_all[, .(
  n_draws_evaluated    = .N,
  mean_marginal_mean   = mean(marginal_mean, na.rm = TRUE),
  median_marginal_mean = median(marginal_mean, na.rm = TRUE),
  sd_marginal_mean     = sd(marginal_mean, na.rm = TRUE),
  ci_2.5               = quantile(marginal_mean, 0.025, na.rm = TRUE),
  ci_97.5              = quantile(marginal_mean, 0.975, na.rm = TRUE),
  prob_positive        = mean(marginal_mean > 0, na.rm = TRUE),
  sign_stable          = (all(marginal_mean > 0) || all(marginal_mean < 0))
), by = .(metric, reg)]

saveRDS(draws_dt_all, file.path(out_causal_dir, "q1_draws_all_estimates.rds"))
saveRDS(q1_robustness, file.path(out_causal_dir, "q1_draws_robustness.rds"))
saveRDS(regional_dt_all, file.path(out_causal_dir, "q1_draws_regional_estimates.rds"))
saveRDS(regional_summary, file.path(out_causal_dir, "q1_draws_regional_summary.rds"))