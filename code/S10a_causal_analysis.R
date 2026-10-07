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

# --- 3. Analysis Part 1: Posterior Mean (Reference Evaluation) ---
cat("\n--- Evaluating Posterior Mean Metrics ---\n")

if(mode == "mean"){
  
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

# --- 4. Analysis Part 2: Posterior Draws (Robustness Check) ---
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
  
# Pre-allocate collector list using total combinations
total_iterations <- length(metrics_list) * length(draw_ids)
draws_results <- vector("list", total_iterations)
counter <- 1
  
for (m in metrics_list) {
  # m <- "E_FRic_loss"
  sub_m <- dt_draws[metrics == m]
  
  for (d in draw_ids) {
    # d <- 922
    sub_md <- sub_m[posterior == d]
    if (nrow(sub_md) == 0) next
      
    fit_d <- tryCatch({
      feols(
        val_mean ~ 1 + i(reg) + i(trees) + i(year),
        data    = sub_md,
        cluster = ~group_number#, 
        #weights = ~weigth_sd
      )
    }, error = function(e) {
      message(sprintf("Model failed for metric %s (draw %s): %s", m, d, e$message))
      return(NULL)
    })
      
    if (!is.null(fit_d)) {
      int_c <- as.numeric(coef(fit_d)["(Intercept)"])
      int_s <- as.numeric(se(fit_d)["(Intercept)"])
      int_p <- as.numeric(pvalue(fit_d)["(Intercept)"])
      ci_d  <- as.numeric(confint(fit_d)["(Intercept)", ])
      
      row_draw <- data.table(
        metric    = as.character(m),
        posterior = as.integer(d),
        n_stands  = as.integer(nrow(sub_md)),
        intercept = int_c,
        se        = int_s,
        t_stat    = int_c / int_s,
        p_value   = int_p,
        ci_2.5    = ci_d[1],
        ci_97.5   = ci_d[2]
      )
      
      draws_results[[counter]] <- row_draw
      
      counter <- counter + 1
    }
  }
}
  
# Remove unused pre-allocated slots and bind
draws_dt_all <- rbindlist(draws_results[seq_len(counter - 1)])
  
# Posterior Credible Summaries across draws (evaluating Q1 robustness)
q1_robustness <- draws_dt_all[, .(
  n_draws_evaluated = .N,
  mean_intercept    = mean(intercept, na.rm = TRUE),
  median_intercept  = median(intercept, na.rm = TRUE),
  sd_intercept      = sd(intercept, na.rm = TRUE),
  ci_2.5            = quantile(intercept, 0.025, na.rm = TRUE),
  ci_97.5           = quantile(intercept, 0.975, na.rm = TRUE),
  prob_positive     = mean(intercept > 0, na.rm = TRUE),
  sign_stable       = (all(intercept > 0) || all(intercept < 0))
), by = .(metric)]
  

saveRDS(draws_dt_all, file.path(out_causal_dir, "q1_draws_all_estimates.rds"))
saveRDS(q1_robustness, file.path(out_causal_dir, "q1_draws_robustness.rds"))
