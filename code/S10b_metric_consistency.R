#' DESCRIPTION:
#' S10b: Evaluation of metric consistency across stand-level predictions
#' restricted to within matched groups (strata).
#' 
#' METHODOLOGY:
#' - Aggregates posterior draws across stand-level metric predictions within groups.
#' - Spearman rank correlation and top 10% overlap evaluated within matched groups.
#' - Posterior probability of opposite directional signals (discordance) per stand.
#' - Group-level summaries and meta-aggregation across groups.
#' 
#' INPUTS:
#' - Group Parquet files: results/metrics/[expected_type]/[run_mode]/metrics_group_*.parquet
#' 
#' OUTPUTS:
#' - Within-Group Rank & Concordance Summaries: results/consistency/q2_within_stratum_*.csv
#' - Stand-Level Discordance Summaries: results/consistency/q2_stand_discordance*.csv
#' - Heatmaps: results/consistency/q2_*_heatmap.png / .pdf

library(tidyverse)
library(arrow)
library(here)
library(patchwork)

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

run_mode <- "draws"

draws_dir <- if (test) {
  file.path(wr_dir, "metrics", "test", expected_string, "draws")
} else {
  file.path(wr_dir, "metrics", expected_string, "draws")
}

output_dir <- file.path(wr_dir, "consistency")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

priority_fraction    <- 0.10
support_threshold    <- 0.90
min_draw_fraction    <- 0.80
min_stands_per_group <- 10L

footprints <- c("E_MSA_loss", "E_PDF_unforced", "E_PDF_forced")
ecological <- c("EG_loss", "E_FRic_loss", "E_RaoQ_loss", "E_BCD")
metrics    <- c(footprints, ecological)

epsilon <- setNames(rep(0, length(metrics)), metrics)

# --- 2. Load Metric Data ---
ds_draws <- arrow::open_dataset(draws_dir)

wide <- ds_draws |>
  select(standid_treated, group_number, posterior, metrics, val_mean) |>
  collect() |>
  pivot_wider(
    names_from = metrics,
    values_from = val_mean
  )

draws <- sort(unique(wide$posterior))
required_draws <- ceiling(min_draw_fraction * length(draws))

all_pairs <- as_tibble(t(combn(metrics, 2L))) |>
  rename(metric_a = V1, metric_b = V2)

comparison_pairs <- tidyr::crossing(metric_a = footprints, metric_b = ecological)
footprint_pairs  <- as_tibble(t(combn(footprints, 2L))) |>
  rename(metric_a = V1, metric_b = V2)

signed_pairs <- bind_rows(
  comparison_pairs |> filter(metric_b != "E_BCD"),
  footprint_pairs
) |> 
  distinct()

safe_q <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  unname(stats::quantile(x, probs = p, na.rm = TRUE))
}

safe_prop <- function(num, den) {
  if (!is.finite(den) || den <= 0L) NA_real_ else num / den
}

get_priority <- function(x, fraction) {
  n <- length(x)
  if (!n) return(logical(0))
  rank(-x, ties.method = "min") <= ceiling(fraction * n)
}

# --- 3. Analysis Part 1: Within-Group Rank Concordance & Priority Overlap ---
cat("\n--- Evaluating Within-Group Rank Concordance & Overlap ---\n")

group_rank_results <- purrr::pmap_dfr(all_pairs, function(metric_a, metric_b) {
  wide |>
    group_by(group_number, draw = posterior) |>
    reframe({
      x  <- .data[[metric_a]]
      y  <- .data[[metric_b]]
      ok <- is.finite(x) & is.finite(y)
      n  <- sum(ok)
      
      if (n < min_stands_per_group || n_distinct(x[ok]) < 2L || n_distinct(y[ok]) < 2L) {
        tibble(
          n_stands       = n,
          rho            = NA_real_,
          n_priority_a   = NA_integer_,
          n_priority_b   = NA_integer_,
          top_overlap    = NA_real_,
          chance_overlap = NA_real_,
          top_lift       = NA_real_
        )
      } else {
        rho_val <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
        pa      <- get_priority(x[ok], priority_fraction)
        pb      <- get_priority(y[ok], priority_fraction)
        ov      <- safe_prop(sum(pa & pb), sum(pa))
        ch      <- mean(pb)
        lift    <- if (is.finite(ch) && ch > 0) ov / ch else NA_real_
        
        tibble(
          n_stands       = n,
          rho            = rho_val,
          n_priority_a   = sum(pa),
          n_priority_b   = sum(pb),
          top_overlap    = ov,
          chance_overlap = ch,
          top_lift       = lift
        )
      }
    }) |>
    mutate(metric_a = metric_a, metric_b = metric_b)
})

readr::write_csv(group_rank_results, file.path(output_dir, "q2_within_stratum_by_draw.csv"))

# Summarize across posterior draws within each group
group_summary <- group_rank_results |>
  group_by(group_number, metric_a, metric_b) |>
  summarise(
    n_treated_stands      = max(n_stands),
    n_usable_draws        = sum(is.finite(rho)),
    group_median_rho      = safe_q(rho, 0.5),
    group_lower_rho_95    = safe_q(rho, 0.025),
    group_upper_rho_95    = safe_q(rho, 0.975),
    group_median_overlap  = safe_q(top_overlap, 0.5),
    group_median_lift     = safe_q(top_lift, 0.5),
    .groups               = "drop"
  )

readr::write_csv(group_summary, file.path(output_dir, "q2_within_stratum_summary.csv"))

# Meta-aggregation across matched groups (without variable masking)
overall_rank_summary <- group_summary |>
  filter(is.finite(group_median_rho)) |>
  group_by(metric_a, metric_b) |>
  summarise(
    n_groups           = n(),
    median_rho         = safe_q(group_median_rho, 0.5),
    lower_rho_95       = safe_q(group_median_rho, 0.025),
    upper_rho_95       = safe_q(group_median_rho, 0.975),
    median_top_overlap = safe_q(group_median_overlap, 0.5),
    lower_top_overlap  = safe_q(group_median_overlap, 0.025),
    upper_top_overlap  = safe_q(group_median_overlap, 0.975),
    median_top_lift    = safe_q(group_median_lift, 0.5),
    .groups            = "drop"
  )

readr::write_csv(overall_rank_summary, file.path(output_dir, "q2_rank_summary.csv"))

# --- 4. Analysis Part 2: Stand-Level Discordance ---
cat("\n--- Evaluating Stand-Level Discordance Probability ---\n")

stand_results <- purrr::pmap_dfr(signed_pairs, function(metric_a, metric_b) {
  wide |>
    group_by(standid_treated, group_number) |>
    summarise({
      val_a      <- .data[[metric_a]]
      val_b      <- .data[[metric_b]]
      eps_a      <- epsilon[[metric_a]]
      eps_b      <- epsilon[[metric_b]]
      
      valid      <- is.finite(val_a) & is.finite(val_b)
      decisive   <- valid & abs(val_a) > eps_a & abs(val_b) > eps_b
      n_opposite <- sum(decisive & (val_a * val_b < 0), na.rm = TRUE)
      
      n_valid    <- sum(valid)
      n_decisive <- sum(decisive)
      
      tibble(
        n_draws_valid             = n_valid,
        n_draws_decisive          = n_decisive,
        n_zero_either             = n_valid - n_decisive,
        fraction_zero_either      = safe_prop(n_valid - n_decisive, n_valid),
        n_opposite                = n_opposite,
        prob_discordance          = safe_prop(n_opposite, n_valid),
        prob_discordance_decisive = safe_prop(n_opposite, n_decisive),
        prob_a_loss               = safe_prop(sum(valid & val_a > eps_a), n_valid),
        prob_b_loss               = safe_prop(sum(valid & val_b > eps_b), n_valid)
      )
    }, .groups = "drop") |>
    mutate(metric_a = metric_a, metric_b = metric_b)
}) |>
  mutate(
    strong_discordance = coalesce(
      n_draws_valid >= required_draws & prob_discordance >= support_threshold,
      FALSE
    ),
    strong_discordance_decisive = coalesce(
      n_draws_decisive >= required_draws & prob_discordance_decisive >= support_threshold,
      FALSE
    )
  ) |>
  relocate(
    standid_treated, group_number, metric_a, metric_b,
    n_draws_valid, n_draws_decisive, n_zero_either, fraction_zero_either,
    n_opposite, prob_discordance, prob_discordance_decisive,
    prob_a_loss, prob_b_loss, strong_discordance, strong_discordance_decisive
  )

stand_summary <- stand_results |>
  group_by(metric_a, metric_b) |>
  summarise(
    n_stands                               = n(),
    n_strong_discordance                   = sum(strong_discordance, na.rm = TRUE),
    proportion_strong_discordance          = mean(strong_discordance, na.rm = TRUE),
    n_strong_discordance_decisive          = sum(strong_discordance_decisive, na.rm = TRUE),
    proportion_strong_discordance_decisive = mean(strong_discordance_decisive, na.rm = TRUE),
    median_prob_discordance                = safe_q(prob_discordance, 0.5),
    median_prob_discordance_decisive       = safe_q(prob_discordance_decisive, 0.5),
    median_fraction_zero_either            = safe_q(fraction_zero_either, 0.5),
    .groups                                = "drop"
  )

readr::write_csv(stand_results, file.path(output_dir, "q2_stand_discordance.csv"))
readr::write_csv(stand_summary, file.path(output_dir, "q2_stand_discordance_summary.csv"))