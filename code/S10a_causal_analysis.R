#' DESCRIPTION:
#' S10a: Evaluation of METSO conservation impact across biodiversity metrics
#' using fixed-effect Ordinary Least Squares (fixest::feols) across W0-W3 weighting schemes.
#' 
#' METHODOLOGY:
#' - Evaluates whether the conditional mean difference (val_mean) between METSO 
#'   and matched BAU controls significantly differs from zero across strata.
#' - Specification: val_mean ~ 1 + i(reg) + i(trees) + i(year), cluster = ~group_number
#' - Weighting Schemes:
#'   - W0 (unweight): unweighted OLS
#'   - W1 (variance): inverse-variance weighted OLS (1 / val_sd^2)
#'   - W2 (domain_variance): variance weighted by MOP continuous distance and out-of-range penalties
#'   - W3 (indomain_variance): inverse-variance weighted OLS restricted to in-domain stands (strict_extrap == 0)
#' - Robustness: Propagates parameter and prediction uncertainty across posterior draws.
#' - Draw-level national and region-standardized marginal means from weighted OLS.
#' 
#' INPUTS:
#' - Group Parquet files: results/metrics/[expected_type]/draws/metrics_group_*.parquet
#' - Group Matches: 
#'  - results/matched_pairs_groups_metso.rds
#'  - results/matched_pairs_groups_control.rds
#'  - results/mop/mop_diagnostics_mask.csv
#' 
#' OUTPUTS:
#' - Summary Tables and Draws in:
#'   - results/causal_analysis/unweight/
#'   - results/causal_analysis/variance/
#'   - results/causal_analysis/domain_variance/
#'   - results/causal_analysis/indomain_variance/
#' - Weight Sensitivity Diagnostics:
#'   - results/causal_analysis/weight_sensitivity_diagnostics.csv

library(tidyverse)
library(arrow)
library(fixest)
library(here)
library(data.table)

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
  stop("Legacy 'mean' mode was removed to avoid metric bias (e.g. FRic). Use 'draws' or 'all' modes only.")
}

draws_dir <- if (test) {
  file.path(wr_dir, "metrics", "test", expected_string, "draws")
} else {
  file.path(wr_dir, "metrics", expected_string, "draws")
}

# --- 2. Load Stratum Metadata & MOP Diagnostics for Weighting ---
matches_metso <- readRDS(file.path(wr_dir, "matched_pairs_groups_metso.rds"))
matches_control <- readRDS(file.path(wr_dir, "matched_pairs_groups_control.rds"))

groups_meta <- bind_rows(matches_metso, matches_control) |>
  select(group_number, reg, trees, year) |>
  distinct() |>
  mutate(
    reg = as.factor(reg),
    trees = as.factor(trees),
    year = as.factor(year)
  )

# Load MOP diagnostics and compute average-approach domain weights
mop_dt <- read_csv(
  file.path(wr_dir, "mop", "mop_diagnostics_mask.csv"),
  show_col_types = FALSE
)

metso_mop <- mop_dt |>
  filter(group == "METSO") |>
  transmute(
    standid_treated = stand_id,
    mop_dist_metso = mop_distance,
    n_out_metso = number_out_vars,
    strict_extrap_metso = mop_strict_extrap
  )

control_mop <- mop_dt |>
  filter(group == "BAU") |>
  left_join(
    matches_control |>
      select(standid, group_number),
    by = c("stand_id" = "standid")
  )

bau_group_mop <- control_mop |>
  group_by(group_number) |>
  summarise(
    mean_bau_dist = mean(mop_distance, na.rm = TRUE),
    mean_bau_out = mean(number_out_vars, na.rm = TRUE),
    .groups = "drop"
  )

metso_group_mop <- matches_metso |>
  transmute(
    standid_treated = standid,
    group_number
  ) |>
  left_join(
    metso_mop,
    by = "standid_treated"
  )

mop_combined <- metso_group_mop |>
  left_join(
    bau_group_mop,
    by = "group_number"
  )

d0 <- 3.0556
d95 <- 5.1517

mop_combined <- mop_combined |>
  mutate(
    bar_D = (mop_dist_metso + mean_bau_dist) / 2,
    bar_K = (n_out_metso + mean_bau_out) / 2
  ) |>
  mutate(
    phi_dist = exp(-pmax(0, bar_D - d0) / (d95 - d0)),
    phi_out = 0.5^bar_K
  )

# --- 3. Load Draws Data & Construct Target Populations ---
cat("\n--- Loading Posterior Draws Data ---\n")

ds_draws <- arrow::open_dataset(draws_dir)

dt_draws <- ds_draws |>
  select(standid_treated, group_number, posterior, metrics, val_mean, val_sd) |>
  collect()

# Floor val_sd == 0 at the 0.1st percentile of positive values to avoid Inf weights
dt_draws <- dt_draws |>
  group_by(metrics) |>
  mutate(
    val_sd_safe = ifelse(val_sd == 0 | is.na(val_sd), quantile(val_sd[val_sd > 0], 0.001, na.rm = TRUE), val_sd),
    weigth_sd = (1 / val_sd_safe)^2
  ) |>
  ungroup()

dt_draws <- dt_draws |>
  left_join(
    groups_meta,
    by = "group_number"
  ) |>
  left_join(
    mop_combined |>
      select(
        standid_treated,
        strict_extrap_metso,
        phi_dist,
        phi_out
      ),
    by = "standid_treated"
  ) |>
  mutate(
    weight_domain = weigth_sd * phi_dist * phi_out
  )

metrics_list <- sort(unique(dt_draws$metrics))
draw_ids <- sort(unique(dt_draws$posterior))

# -------------------------------------------------------------------------
# Target Population Grids:
# 1) National Grid: Represents the full observed METSO population (all 29,159 stands).
#    Weights strata by their actual share in the entire national conservation portfolio.
# 2) In-Domain Grid: Represents only stands strictly within the FBS environmental envelope.
#    Reweights the target population to the stratum composition of non-extrapolated stands.
# -------------------------------------------------------------------------

# Full National Population Grid
stand_pop <- dt_draws |>
  select(standid_treated, group_number, reg, trees, year) |>
  distinct()

national_grid <- stand_pop |>
  count(reg, trees, year, name = "n_target") |>
  mutate(
    share = n_target / sum(n_target)
  )

tree_year_mix <- stand_pop |>
  count(trees, year, name = "n_target") |>
  mutate(
    share = n_target / sum(n_target)
  )

regional_grid <- tidyr::crossing(
  reg = unique(stand_pop$reg),
  tree_year_mix
)

regional_counts <- stand_pop |>
  group_by(reg) |>
  summarise(
    n_region_stands = n(),
    n_region_groups = n_distinct(group_number),
    .groups = "drop"
  )

# In-Domain Population Grid (W3)
stand_pop_indomain <- dt_draws |>
  filter(strict_extrap_metso == 0) |>
  select(standid_treated, group_number, reg, trees, year) |>
  distinct()

national_grid_indomain <- stand_pop_indomain |>
  count(reg, trees, year, name = "n_target") |>
  mutate(
    share = n_target / sum(n_target)
  )

tree_year_mix_indomain <- stand_pop_indomain |>
  count(trees, year, name = "n_target") |>
  mutate(
    share = n_target / sum(n_target)
  )

regional_grid_indomain <- tidyr::crossing(
  reg = unique(stand_pop_indomain$reg),
  tree_year_mix_indomain
)

regional_counts_indomain <- stand_pop_indomain |>
  group_by(reg) |>
  summarise(
    n_region_stands = n(),
    n_region_groups = n_distinct(group_number),
    .groups = "drop"
  )

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
  l <- colSums(x * w)
  estimate_from_coef <- sum(l * b)
  if (!isTRUE(all.equal(estimate, estimate_from_coef, tolerance = 1e-10))) {
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

# --- 4. Loop over Weighting Schemes W0 to W3 ---
weight_schemes <- c(
  "unweight",
  "variance",
  "domain_variance",
  "indomain_variance"
)

for (ws in weight_schemes) {
  cat(sprintf("\n--- Evaluating Weighting Scheme: %s ---\n", ws))
  
  out_causal_dir <- file.path(wr_dir, "causal_analysis", ws)
  dir.create(out_causal_dir, recursive = TRUE, showWarnings = FALSE)
  
  total_iterations <- length(metrics_list) * length(draw_ids)
  draws_results <- vector("list", total_iterations)
  regional_results <- vector("list", total_iterations)
  counter <- 1
  
  for (m in metrics_list) {
    sub_m <- dt_draws |>
      filter(metrics == m)
    
    for (d in draw_ids) {
      sub_md <- sub_m |>
        filter(posterior == d)
      
      if (nrow(sub_md) == 0) next
      
      # Select population reference grid based on the evaluated domain
      if (ws == "indomain_variance") {
        sub_md <- sub_md |>
          filter(strict_extrap_metso == 0)
        national_grid_current <- national_grid_indomain
        regional_grid_current <- regional_grid_indomain
        regional_counts_current <- regional_counts_indomain
        current_pop <- stand_pop_indomain
        
        if (nrow(sub_md) == 0) next
      } else {
        national_grid_current <- national_grid
        regional_grid_current <- regional_grid
        regional_counts_current <- regional_counts
        current_pop <- stand_pop
      }
      
      fit_d <- tryCatch({
        if (ws == "unweight") {
          feols(
            val_mean ~ 1 + i(reg) + i(trees) + i(year),
            data = sub_md,
            cluster = ~group_number
          )
        } else if (ws %in% c("variance", "indomain_variance")) {
          feols(
            val_mean ~ 1 + i(reg) + i(trees) + i(year),
            data = sub_md,
            cluster = ~group_number,
            weights = ~weigth_sd
          )
        } else if (ws == "domain_variance") {
          feols(
            val_mean ~ 1 + i(reg) + i(trees) + i(year),
            data = sub_md,
            cluster = ~group_number,
            weights = ~weight_domain
          )
        }
      }, error = function(e) {
        message(sprintf(
          "Model failed for metric %s (draw %s, scheme %s): %s",
          m,
          d,
          ws,
          e$message
        ))
        return(NULL)
      })
      
      if (!is.null(fit_d)) {
        national_margin <- marginal_estimate(
          fit_d,
          national_grid_current
        )
        
        row_draw <- bind_cols(
          tibble(
            metric = as.character(m),
            posterior = as.integer(d),
            n_stands = nrow(current_pop),
            n_groups = n_distinct(sub_md$group_number),
            n_fitted = nobs(fit_d)
          ),
          national_margin
        )
        
        region_margin <- regional_grid_current |>
          group_split(reg, .keep = TRUE) |>
          map_dfr(function(g) {
            bind_cols(
              tibble(reg = g$reg[1]),
              marginal_estimate(fit = fit_d, grid = g)
            )
          })
        
        region_margin <- region_margin |>
          left_join(
            regional_counts_current,
            by = "reg"
          ) |>
          rename(
            marginal_mean = estimate
          ) |>
          mutate(
            metric = as.character(m),
            posterior = as.integer(d),
            t_stat = marginal_mean / se,
            p_value = 2 * pnorm(abs(t_stat), lower.tail = FALSE)
          ) |>
          select(
            metric,
            posterior,
            reg,
            n_region_stands,
            n_region_groups,
            marginal_mean,
            se,
            t_stat,
            p_value,
            ci_2.5,
            ci_97.5
          )
        
        draws_results[[counter]] <- row_draw
        regional_results[[counter]] <- region_margin
        
        counter <- counter + 1
      }
    }
  }
  
  draws_dt_all <- draws_results |>
    keep(Negate(is.null)) |>
    bind_rows()
  
  regional_dt_all <- regional_results |>
    keep(Negate(is.null)) |>
    bind_rows()
  
  q1_robustness <- draws_dt_all |>
    group_by(metric) |>
    summarise(
      n_draws_evaluated = n(),
      estimate_mean = mean(estimate, na.rm = TRUE),
      estimate_median = median(estimate, na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      ci_2.5 = quantile(estimate, 0.025, na.rm = TRUE),
      ci_97.5 = quantile(estimate, 0.975, na.rm = TRUE),
      prob_positive = mean(estimate > 0, na.rm = TRUE),
      sign_stable = all(estimate > 0) || all(estimate < 0),
      .groups = "drop"
    )
  
  regional_summary <- regional_dt_all |>
    group_by(metric, reg) |>
    summarise(
      n_draws_evaluated = n(),
      mean_marginal_mean = mean(marginal_mean, na.rm = TRUE),
      median_marginal_mean = median(marginal_mean, na.rm = TRUE),
      sd_marginal_mean = sd(marginal_mean, na.rm = TRUE),
      ci_2.5 = quantile(marginal_mean, 0.025, na.rm = TRUE),
      ci_97.5 = quantile(marginal_mean, 0.975, na.rm = TRUE),
      prob_positive = mean(marginal_mean > 0, na.rm = TRUE),
      sign_stable = all(marginal_mean > 0) || all(marginal_mean < 0),
      .groups = "drop"
    )
  
  write_csv(
    draws_dt_all,
    file.path(out_causal_dir, "q1_draws_all_estimates.csv")
  )
  write_csv(
    q1_robustness,
    file.path(out_causal_dir, "q1_draws_robustness.csv")
  )
  write_csv(
    regional_dt_all,
    file.path(out_causal_dir, "q1_draws_regional_estimates.csv")
  )
  write_csv(
    regional_summary,
    file.path(out_causal_dir, "q1_draws_regional_summary.csv")
  )
  
  cat(sprintf("Completed weighting scheme: %s\n", ws))
}

# --- 5. Weight Diagnostics and Sensitivity Checks ---
cat("\n--- Running Weight Sensitivity and Finite Value Diagnostics ---\n")

# Distribution of weights and MOP reliability penalty factors
weight_distribution_summary <- dt_draws |>
  summarise(
    n = n(),
    min_weigth_sd = min(weigth_sd, na.rm = TRUE),
    q25_weigth_sd = quantile(weigth_sd, 0.25, na.rm = TRUE),
    median_weigth_sd = median(weigth_sd, na.rm = TRUE),
    q75_weigth_sd = quantile(weigth_sd, 0.75, na.rm = TRUE),
    max_weigth_sd = max(weigth_sd, na.rm = TRUE),
    min_phi_dist = min(phi_dist, na.rm = TRUE),
    median_phi_dist = median(phi_dist, na.rm = TRUE),
    max_phi_dist = max(phi_dist, na.rm = TRUE),
    min_phi_out = min(phi_out, na.rm = TRUE),
    median_phi_out = median(phi_out, na.rm = TRUE),
    max_phi_out = max(phi_out, na.rm = TRUE),
    min_weight_domain = min(weight_domain, na.rm = TRUE),
    median_weight_domain = median(weight_domain, na.rm = TRUE),
    max_weight_domain = max(weight_domain, na.rm = TRUE)
  )

# Verify zero counts of non-finite or infinite weights post-flooring
weight_integrity_check <- dt_draws |>
  summarise(
    n_total = n(),
    n_zero_sd = sum(val_sd == 0, na.rm = TRUE),
    n_missing_sd = sum(is.na(val_sd)),
    n_nonfinite_sd = sum(!is.finite(val_sd)),
    n_infinite_weigth_sd = sum(is.infinite(weigth_sd)),
    n_nonfinite_weigth_sd = sum(!is.finite(weigth_sd)),
    n_infinite_weight_domain = sum(is.infinite(weight_domain)),
    n_nonfinite_weight_domain = sum(!is.finite(weight_domain))
  )

# Evaluate weight dispersion (CV) and correlation between weights and effect sizes across posterior draws
weight_variation_by_draw <- dt_draws |>
  group_by(metrics, posterior) |>
  summarise(
    n = n(),
    cv_weigth_sd = sd(weigth_sd, na.rm = TRUE) / mean(weigth_sd, na.rm = TRUE),
    cv_weight_domain = sd(weight_domain, na.rm = TRUE) / mean(weight_domain, na.rm = TRUE),
    cor_weight_outcome = cor(weight_domain, val_mean, use = "complete.obs"),
    .groups = "drop"
  )

weight_sensitivity_summary <- weight_variation_by_draw |>
  group_by(metrics) |>
  summarise(
    median_cv_weigth_sd = median(cv_weigth_sd, na.rm = TRUE),
    median_cv_weight_domain = median(cv_weight_domain, na.rm = TRUE),
    median_cor_weight_outcome = median(cor_weight_outcome, na.rm = TRUE),
    min_cor_weight_outcome = min(cor_weight_outcome, na.rm = TRUE),
    max_cor_weight_outcome = max(cor_weight_outcome, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(abs(median_cor_weight_outcome)))

print(weight_distribution_summary)
print(weight_integrity_check)
print(weight_sensitivity_summary)

write_csv(
  weight_sensitivity_summary,
  file.path(wr_dir, "causal_analysis", "weight_sensitivity_diagnostics.csv")
)

cat("Diagnostics successfully executed and saved.\n")