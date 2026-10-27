#' DESCRIPTION:
#' S10a: Evaluation of METSO conservation impact across biodiversity metrics
#' using fixed-effect Ordinary Least Squares (fixest::feols) across W0-W3 weighting schemes.
#' 
#' METHODOLOGY:
#' - Evaluates whether the conditional mean difference (val_mean) between METSO 
#'   and matched BAU controls significantly differs from zero across groups.
#' - Specification: val_mean ~ 1 + i(reg) + i(trees) + i(year), cluster = ~group_number
#' - Weighting Schemes:
#'   - W0 (unweight): unweighted OLS across all eligible METSO stands.
#'   - W1 (variance): inverse-variance weighted OLS (1 / val_sd^2) across all eligible METSO stands.
#'   - W2 (domain_variance): variance weighted by continuous MOP distance and extrapolation penalties.
#'   - W3 (indomain_variance): inverse-variance weighted OLS restricted strictly to in-domain METSO
#'     stands evaluated against in-domain BAU controls (from S09 comparison_mode == "indomain_both").
#' - Pre-estimation exclusion: Logs non-finite, missing, or zero-variance weights to CSV before fitting.
#' - Robustness: Propagates parameter and prediction uncertainty across posterior draws.
#' - Draw-level national and region-standardized marginal means from weighted OLS.
#' 
#' INPUTS:
#' - Full Group Parquet files: results/metrics/[expected_type]/draws/metrics_group_*.parquet
#' - In-Domain Group Parquet files: results/metrics_indomain_both/[expected_type]/draws/metrics_group_*.parquet
#' - Group Matches: 
#'   - results/matched_pairs_groups_metso.rds
#'   - results/matched_pairs_groups_control.rds
#'   - results/mop/mop_diagnostics_mask.csv
#' 
#' OUTPUTS:
#' - Summary Tables and Draws in:
#'   - results/causal_analysis/unweight/
#'   - results/causal_analysis/variance/
#'   - results/causal_analysis/domain_variance/
#'   - results/causal_analysis/indomain_variance/
#' - Exclusion and Sensitivity Logs:
#'   - results/causal_analysis/weight_exclusion_log.csv
#'   - results/causal_analysis/weight_sensitivity_diagnostics.csv

library(tidyverse)
library(arrow)
library(fixest)
library(here)
library(data.table)

source(file.path("code", "config_model.R"))

message("Start time: ", Sys.time())

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

run_mode <- "draws"
if (run_mode == "mean") {
  stop("Legacy 'mean' mode was removed to avoid metric bias (e.g. FRic). Use 'draws' or 'all' modes only.")
}

draws_dir_full <- if (test) {
  file.path(wr_dir, "metrics", "test", expected_string, "draws")
} else {
  file.path(wr_dir, "metrics", expected_string, "draws")
}

draws_dir_indomain <- if (test) {
  file.path(wr_dir, "metrics_indomain_both", "test", expected_string, "draws")
} else {
  file.path(wr_dir, "metrics_indomain_both", expected_string, "draws")
}

out_causal_base <- file.path(wr_dir, "causal_analysis")
dir.create(out_causal_base, recursive = TRUE, showWarnings = FALSE)


# --- 2. Load groups Metadata & MOP Continuous Diagnostics ---
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
    matches_control |> select(standid, group_number),
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
  left_join(metso_mop, by = "standid_treated")

mop_combined <- metso_group_mop |>
  left_join(bau_group_mop, by = "group_number")

# Continuous MOP penalties for W2
d0 <- median(mop_dt$mop_distance, na.rm = TRUE)
d95 <- quantile(mop_dt$mop_distance, 0.95, na.rm = TRUE)

mop_combined <- mop_combined |>
  mutate(
    bar_D = (mop_dist_metso + mean_bau_dist) / 2,
    bar_K = (n_out_metso + mean_bau_out) / 2,
    phi_dist = exp(-pmax(0, bar_D - d0) / (d95 - d0)),
    phi_out = 0.5^bar_K
  )


# --- 3. Load Datasets & Build Population Grids ---
message("Loading full and in-domain metric datasets...")

ds_full <- arrow::open_dataset(draws_dir_full)
dt_full <- ds_full |>
  select(standid_treated, group_number, posterior, metrics, val_mean, val_sd, n_controls) |>
  collect() |>
  left_join(groups_meta, by = "group_number") |>
  left_join(
    mop_combined |> select(standid_treated, strict_extrap_metso, phi_dist, phi_out),
    by = "standid_treated"
  )

ds_indomain <- arrow::open_dataset(draws_dir_indomain)
dt_indomain <- ds_indomain |>
  select(standid_treated, group_number, posterior, metrics, val_mean, val_sd, n_controls) |>
  collect() |>
  left_join(groups_meta, by = "group_number")


# --- 4. Pre-estimation Weight Calculation & Stand Exclusion ---
# Compute inverse-variance weights
dt_full <- dt_full |>
  mutate(
    raw_var = val_sd^2,
    weight_sd = 1 / raw_var,
    weight_domain = weight_sd * phi_dist * phi_out
  )

dt_indomain <- dt_indomain |>
  mutate(
    raw_var = val_sd^2,
    weight_sd = 1 / raw_var
  )

# Flag and log excluded records
exclusions_full <- dt_full |>
  filter(is.na(val_mean) | is.na(weight_sd) | !is.finite(weight_sd) | weight_sd <= 0) |>
  transmute(
    dataset = "full",
    metrics,
    posterior,
    standid_treated,
    group_number,
    reg,
    trees,
    year,
    val_mean,
    val_sd,
    n_controls,
    reason = case_when(
      is.na(val_mean) ~ "missing_val_mean",
      is.na(val_sd) ~ "missing_val_sd",
      val_sd == 0 ~ "zero_val_sd_infinite_weight",
      !is.finite(weight_sd) ~ "non_finite_weight",
      weight_sd <= 0 ~ "non_positive_weight",
      TRUE ~ "other"
    )
  )

exclusions_indomain <- dt_indomain |>
  filter(is.na(val_mean) | is.na(weight_sd) | !is.finite(weight_sd) | weight_sd <= 0) |>
  transmute(
    dataset = "indomain_both",
    metrics,
    posterior,
    standid_treated,
    group_number,
    reg,
    trees,
    year,
    val_mean,
    val_sd,
    n_controls,
    reason = case_when(
      is.na(val_mean) ~ "missing_val_mean",
      is.na(val_sd) ~ "missing_val_sd",
      val_sd == 0 ~ "zero_val_sd_infinite_weight",
      !is.finite(weight_sd) ~ "non_finite_weight",
      weight_sd <= 0 ~ "non_positive_weight",
      TRUE ~ "other"
    )
  )

exclusion_log <- bind_rows(exclusions_full, exclusions_indomain)
write_csv(exclusion_log, file.path(out_causal_base, "weight_exclusion_log.csv"))

# Retain strictly finite, valid records for estimation
dt_full_clean <- dt_full |>
  filter(!is.na(val_mean))

dt_indomain_clean <- dt_indomain |>
  filter(!is.na(val_mean) & is.finite(weight_sd) & weight_sd > 0)

# Build Target Population Grids
# National & Regional Grids for Full Cohort (W0, W1, W2)
stand_pop_full <- dt_full_clean |>
  select(standid_treated, group_number, reg, trees, year) |>
  distinct()

national_grid_full <- stand_pop_full |>
  count(reg, trees, year, name = "n_target") |>
  mutate(share = n_target / sum(n_target))

tree_year_mix_full <- stand_pop_full |>
  count(trees, year, name = "n_target") |>
  mutate(share = n_target / sum(n_target))

regional_grid_full <- tidyr::crossing(
  reg = unique(stand_pop_full$reg),
  tree_year_mix_full
)

regional_counts_full <- stand_pop_full |>
  group_by(reg) |>
  summarise(
    n_region_stands = n(),
    n_region_groups = n_distinct(group_number),
    .groups = "drop"
  )

# National & Regional Grids for Restricted In-Domain Cohort (W3)
stand_pop_indomain <- dt_indomain_clean |>
  select(standid_treated, group_number, reg, trees, year) |>
  distinct()

national_grid_indomain <- stand_pop_indomain |>
  count(reg, trees, year, name = "n_target") |>
  mutate(share = n_target / sum(n_target))

tree_year_mix_indomain <- stand_pop_indomain |>
  count(trees, year, name = "n_target") |>
  mutate(share = n_target / sum(n_target))

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


# --- 5. Helper Function ---
marginal_estimate <- function(fit, grid, weight_col = "share") {
  w <- grid[[weight_col]]
  w <- w / sum(w)
  pred <- predict(fit, newdata = as.data.frame(grid))
  estimate <- sum(w * pred)
  b <- coef(fit)
  V <- vcov(fit)
  
  x <- matrix(0, nrow = nrow(grid), ncol = length(b), dimnames = list(NULL, names(b)))
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


# --- 6. Model Estimation Loop Across Weighting Schemes (W0 - W3) ---
weight_schemes <- c(
  "unweight",          # W0
  "variance",          # W1
  "domain_variance",   # W2
  "indomain_variance"  # W3
)

metrics_list <- sort(unique(dt_full$metrics))
draw_ids <- sort(unique(dt_full$posterior))

for (ws in weight_schemes) {
  message("Evaluating Weighting Scheme: ", ws)
  # ws <- "variance"
  out_causal_dir <- file.path(out_causal_base, ws)
  dir.create(out_causal_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Select input dataset, grids, and weights according to scheme
  if (ws == "indomain_variance") {
    dt_active <- dt_indomain_clean
    national_grid_current <- national_grid_indomain
    regional_grid_current <- regional_grid_indomain
    regional_counts_current <- regional_counts_indomain
    current_pop <- stand_pop_indomain
  } else {
    dt_active <- dt_full_clean
    if (ws %in% c("variance", "domain_variance")) {
      dt_active <- dt_active |> filter(is.finite(weight_sd) & weight_sd > 0)
    }
    national_grid_current <- national_grid_full
    regional_grid_current <- regional_grid_full
    regional_counts_current <- regional_counts_full
    current_pop <- stand_pop_full
  }
  
  total_iterations <- length(metrics_list) * length(draw_ids)
  draws_results <- vector("list", total_iterations)
  regional_results <- vector("list", total_iterations)
  counter <- 1
  
  for (m in metrics_list) {
    # m <- "E_MSA_loss"
    sub_m <- dt_active |> filter(metrics == m)
    
    for (d in draw_ids) {
      # d <- 164
      sub_md <- sub_m |> filter(posterior == d)
      if (nrow(sub_md) == 0) next
      
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
            weights = ~weight_sd
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
        message(sprintf("Model failed for metric %s (draw %s, scheme %s): %s", m, d, ws, e$message))
        return(NULL)
      })
      
      if (!is.null(fit_d)) {
        national_margin <- marginal_estimate(fit = fit_d, grid = national_grid_current)
        
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
          }) |>
          left_join(regional_counts_current, by = "reg") |>
          rename(marginal_mean = estimate) |>
          mutate(
            metric = as.character(m),
            posterior = as.integer(d),
            t_stat = marginal_mean / se,
            p_value = 2 * pnorm(abs(t_stat), lower.tail = FALSE)
          ) |>
          select(
            metric, posterior, reg, n_region_stands, n_region_groups,
            marginal_mean, se, t_stat, p_value, ci_2.5, ci_97.5
          )
        
        draws_results[[counter]] <- row_draw
        regional_results[[counter]] <- region_margin
        counter <- counter + 1
      }
    }
  }
  
  draws_dt_all <- draws_results |> keep(Negate(is.null)) |> bind_rows()
  regional_dt_all <- regional_results |> keep(Negate(is.null)) |> bind_rows()
  
  # Summarize across posterior draws
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
  
  write_csv(draws_dt_all, file.path(out_causal_dir, "q1_draws_all_estimates.csv"))
  write_csv(q1_robustness, file.path(out_causal_dir, "q1_draws_robustness.csv"))
  write_csv(regional_dt_all, file.path(out_causal_dir, "q1_draws_regional_estimates.csv"))
  write_csv(regional_summary, file.path(out_causal_dir, "q1_draws_regional_summary.csv"))
  
  message("Completed weighting scheme: ", ws)
}


# --- 7. Weight Sensitivity Diagnostics ---
message("Running weight sensitivity diagnostics...")

weight_variation_by_draw <- dt_full_clean |>
  filter(is.finite(weight_sd) & weight_sd > 0) |>
  group_by(metrics, posterior) |>
  summarise(
    n = n(),
    cv_weight_sd = sd(weight_sd, na.rm = TRUE) / mean(weight_sd, na.rm = TRUE),
    cv_weight_domain = sd(weight_domain, na.rm = TRUE) / mean(weight_domain, na.rm = TRUE),
    cor_weight_outcome = cor(weight_domain, val_mean, use = "complete.obs"),
    .groups = "drop"
  )

weight_sensitivity_summary <- weight_variation_by_draw |>
  group_by(metrics) |>
  summarise(
    median_cv_weight_sd = median(cv_weight_sd, na.rm = TRUE),
    median_cv_weight_domain = median(cv_weight_domain, na.rm = TRUE),
    median_cor_weight_outcome = median(cor_weight_outcome, na.rm = TRUE),
    min_cor_weight_outcome = min(cor_weight_outcome, na.rm = TRUE),
    max_cor_weight_outcome = max(cor_weight_outcome, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(abs(median_cor_weight_outcome)))

write_csv(
  weight_sensitivity_summary,
  file.path(out_causal_base, "weight_sensitivity_diagnostics.csv")
)

# --- 8. Consolidate Outputs Across Weighting Schemes ---
message("Consolidating draw-level and summary tables across schemes...")

scheme_mapping <- tibble(
  weighting_scheme = weight_schemes,
  scheme_code = c("W0", "W1", "W2", "W3")
)

# 1. Consolidate draw-level national estimates (q1_draws_all_estimates.csv)
all_draws_national <- map2_dfr(
  scheme_mapping$weighting_scheme,
  scheme_mapping$scheme_code,
  function(ws, code) {
    file_path <- file.path(out_causal_base, ws, "q1_draws_all_estimates.csv")
    read_csv(file_path, show_col_types = FALSE) |>
      mutate(
        weighting_scheme = ws,
        scheme_code = code,
        .before = 1
      )
  }
)
draws_national_out <- file.path(out_causal_base, "all_schemes_national_draws.csv")
write_csv(all_draws_national, draws_national_out)

# 2. Consolidate national robustness summaries (q1_draws_robustness.csv)
all_robustness <- map2_dfr(
  scheme_mapping$weighting_scheme,
  scheme_mapping$scheme_code,
  function(ws, code) {
    file_path <- file.path(out_causal_base, ws, "q1_draws_robustness.csv")
    read_csv(file_path, show_col_types = FALSE) |>
      mutate(
        weighting_scheme = ws,
        scheme_code = code,
        .before = 1
      )
  }
)
robustness_out <- file.path(out_causal_base, "all_schemes_national_robustness_summary.csv")
write_csv(all_robustness, robustness_out)

# 3. Consolidate draw-level regional estimates (q1_draws_regional_estimates.csv)
all_draws_regional <- map2_dfr(
  scheme_mapping$weighting_scheme,
  scheme_mapping$scheme_code,
  function(ws, code) {
    file_path <- file.path(out_causal_base, ws, "q1_draws_regional_estimates.csv")
    read_csv(file_path, show_col_types = FALSE) |>
      mutate(
        weighting_scheme = ws,
        scheme_code = code,
        .before = 1
      )
  }
)
draws_regional_out <- file.path(out_causal_base, "all_schemes_regional_draws.csv")
write_csv(all_draws_regional, draws_regional_out)

# 4. Consolidate regional summaries (q1_draws_regional_summary.csv)
all_regional <- map2_dfr(
  scheme_mapping$weighting_scheme,
  scheme_mapping$scheme_code,
  function(ws, code) {
    file_path <- file.path(out_causal_base, ws, "q1_draws_regional_summary.csv")
    read_csv(file_path, show_col_types = FALSE) |>
      mutate(
        weighting_scheme = ws,
        scheme_code = code,
        .before = 1
      )
  }
)
regional_out <- file.path(out_causal_base, "all_schemes_regional_summary.csv")
write_csv(all_regional, regional_out)
