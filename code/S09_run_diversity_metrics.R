#' DESCRIPTION:
#' S09: Calculates site-level biodiversity metrics and performs fixed-effect
#' regression (fixest::feols) across METSO stands and matched BAU/control groups.
#' 
#' METHODOLOGY:
#' - Modes: Posterior mean ("mean"), n draws ("draws"), or all draws ("all").
#' - Comparison Modes:
#'   - "full": All matched METSO stands against all available matched BAU controls.
#'   - "indomain_both": Restricted strictly to in-domain METSO stands evaluated
#'     against in-domain BAU controls (mop_strict_extrap == 0).
#' - Aggregation: Averages pairwise comparisons against group controls in-memory.
#' - Parallelism: Distributes METSO stand chunks per group across CPU cores.
#' - Statistical Modeling: Evaluates treatment elasticity via fixed-effects OLS.
#' 
#' INPUTS:
#' - Model Metadata: models/unfitted_RData/unfitted_[model_id].RData
#' - Traits Data: data/traits/TrData_complete.rds
#' - Group Matches: results/matched_pairs_groups_metso.rds & _control.rds
#' - MOP Diagnostics: results/mop/mop_diagnostics_mask.csv (when comparison_mode == "indomain_both")
#' - Predictions: results/predictions/ & results/predictions_PA/ (Arrow partition)
#' 
#' OUTPUTS:
#' - Group Parquet files: 
#'   - results/metrics/[expected_type]/metrics_group_[N].parquet (full)
#'   - results/metrics_indomain_both/[expected_type]/metrics_group_[N].parquet (indomain_both)

library(tidyverse)
library(Hmsc)
library(parallel)
library(here)
library(arrow)
library(data.table)
library(FD)

source(file.path("code", "_utilities_diversity_metrics_functions.R"))
source(file.path("code", "config_model.R"))
name_model <- generate_run_name(run_config)

message("    Start time: ", Sys.time())

# --- Configuration ---
# Modes: "draws", "all". Legacy "mean" disabled.
run_mode <- "draws"
n_draws  <- 100
seed     <- 11072024
set.seed(seed)

# Comparison mode: "full" or "indomain_both" (The former for W0-W3 and the latter for W4 restricted pool)
comparison_mode <- "indomain_both"

# Prevent running legacy 'mean' mode for production
if (run_mode == "mean") {
  stop("Legacy 'mean' mode is disabled for production to avoid metric bias (e.g. FRic). Use 'draws' or 'all'.")
}

selected_draws <- NULL
if (run_mode == "draws") {
  selected_draws <- sort(sample.int(4000, n_draws))
}

# Execution mode flag: set to FALSE for full production run
test <- FALSE

# Expected abundance predictions vs stochastic realizations
expected_string <- "_expected_true"

os <- Sys.info()['sysname']
if (os == "Windows") {
  n_cores <- max(1, floor(parallel::detectCores() / 2))
  wr_dir <- "results"
} else if (os %in% c("Linux", "Darwin")) {
  n_cores <- min(32, floor(parallel::detectCores() / 2))
  node <- Sys.info()["nodename"]
  
  if (node == "lema") {
    wr_dir <- NULL
    worker_tmp <- NULL
  } else if (node == "gurnah") {
    wr_dir <- "/export/scratch/tmp/munozcs/results"
  } else {
    wr_dir <- NULL
  }
}

# Resolve directory destination based on comparison mode
metrics_subfolder <- if (comparison_mode == "indomain_both") {
  "metrics_indomain_both"
} else {
  "metrics"
}

output_dir <- if (test) {
  file.path(wr_dir, metrics_subfolder, "test", expected_string, run_mode)
} else {
  file.path(wr_dir, metrics_subfolder, expected_string, run_mode)
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# 1. Load Species and Traits
load(file.path("models", "unfitted_RData", paste0("unfitted_", run_config$model_id, ".RData")))
name_spp <- colnames(models$fbs_M016$Y)

tr_proc_file <- file.path(wr_dir, "metrics", "TrData_processed.rds")
if (!file.exists(tr_proc_file)) {
  cat("----- Calculating TrData processing...\n")
  TrData <- readRDS(file.path("data", "traits", "TrData_complete.rds"))
  rownames(TrData) <- name_spp
  TrData_numeric <- TrData |> 
    select(where(is.numeric), -any_of(c("Total.individuals", "Female", "Male", "Unknown", "Complete.measures")))
  
  dummy_a <- matrix(1, nrow = 1, ncol = length(name_spp))
  colnames(dummy_a) <- name_spp
  
  TrData_processed <- dbFD_preprocess_traits(x = TrData_numeric, a = dummy_a)
  saveRDS(TrData_processed, tr_proc_file)
} else {
  TrData_processed <- readRDS(tr_proc_file)
}

cum_var <- cumsum(TrData_processed$eigenvalues) / sum(TrData_processed$eigenvalues) * 100
TrData_processed$traits.FRic <- TrData_processed$traits.FRic[, cum_var <= 95]
divmax <- ade4::divcmax(as.dist(TrData_processed$x.dist))$value


# 2. Prepare Match and Arrow Connection
matches_metso <- readRDS(file.path(wr_dir, "matched_pairs_groups_metso.rds"))
matches_control <- readRDS(file.path(wr_dir, "matched_pairs_groups_control.rds"))
groups_df <- bind_rows(matches_metso, matches_control)

# Filter stands and evaluate group eligibility according to comparison_mode
if (comparison_mode == "indomain_both") {
  mop_mask_file <- file.path(wr_dir, "mop", "mop_diagnostics_mask.csv")
  if (!file.exists(mop_mask_file)) {
    stop("MOP diagnostics mask file not found: ", mop_mask_file)
  }
  
  mop_dt <- readr::read_csv(mop_mask_file, show_col_types = FALSE)
  
  # Identify strictly in-domain stand IDs
  indomain_stands <- mop_dt |>
    filter(mop_strict_extrap == 0) |>
    pull(stand_id) |>
    as.character()
  
  # Filter both METSO and BAU control stands by the in-domain condition
  groups_df <- groups_df |>
    filter(as.character(standid) %in% indomain_stands)
  
  # Require at least one METSO stand and at least two BAU controls to retain stratum
  group_counts <- groups_df |>
    group_by(group_number, metso) |>
    summarise(n = n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = metso, values_from = n, names_prefix = "m_") |>
    filter(!is.na(m_1) & m_1 >= 1, !is.na(m_0) & m_0 >= 2)
  
  common_groups <- group_counts$group_number
  cat(sprintf("[%s] Mode indomain_both: %d eligible groups retained with >= 2 controls.\n", 
              Sys.time(), length(common_groups)))
} else {
  # Mode "full": Standard intersection of groups with both METSO and control present
  common_groups <- intersect(
    unique(groups_df$group_number[groups_df$metso == 1]),
    unique(groups_df$group_number[groups_df$metso == 0])
  )
  cat(sprintf("[%s] Mode full: %d groups identified.\n", Sys.time(), length(common_groups)))
}

if (test) {
  test_group_file <- file.path(wr_dir, "test_group_stands.rds")
  
  if (file.exists(test_group_file)) {
    cat("Loading pre-defined test strata from test_group_stands.rds...\n")
    test_meta <- readRDS(test_group_file)
    common_groups <- intersect(common_groups, unique(test_meta$group_number))
  } else {
    cat("Selecting small strata groups directly from groups_df...\n")
    max_stands_per_scenario <- 50
    small_groups <- groups_df |>
      filter(group_number %in% common_groups) |>
      group_by(group_number, metso) |>
      summarise(n_stands = n(), .groups = "drop") |>
      pivot_wider(names_from = metso, values_from = n_stands, names_prefix = "m_") |>
      filter(
        m_1 >= 5 & m_1 <= max_stands_per_scenario,
        m_0 >= 5 & m_0 <= max_stands_per_scenario
      ) |>
      pull(group_number)
    
    common_groups <- head(small_groups, 10)
  }
  
  cat(sprintf("Test: processing %d lightweight groups.\n", length(common_groups)))
}

# Hurdle predictions
predY_ds <- open_dataset(
  file.path(wr_dir, "predictions"), 
  hive_style = TRUE, 
  unify_schemas = TRUE
) |> 
  filter(
    modelid == run_config$model_id, 
    expected_type == expected_string
  )

# Parallel PA probability predictions
predPA_ds <- open_dataset(
  file.path(wr_dir, "predictions_PA"), 
  hive_style = TRUE, 
  unify_schemas = TRUE
) |> 
  filter(
    modelid == run_config$model_id, 
    expected_type == expected_string
  )


# 3. Helper Functions
prepare_stand_matrix <- function(sub_df, mode, n_target, cols_drop, selected_draws = NULL) {
  sub_df <- sub_df[order(sub_df$posterior), ]
  mat <- as.matrix(sub_df[, !(names(sub_df) %in% cols_drop)])
  
  if (mode == "mean") {
    mat_eval <- matrix(colMeans(mat, na.rm = TRUE), nrow = 1)
    colnames(mat_eval) <- colnames(mat)
    return(list(mat = mat_eval, draws = "mean"))
  } else if (mode == "draws") {
    if (is.null(selected_draws)) {
      avail <- unique(sub_df$posterior)
      selected_draws <- sort(sample(avail, min(n_target, length(avail))))
    }
    mat_eval <- mat[sub_df$posterior %in% selected_draws, , drop = FALSE]
    return(list(mat = mat_eval, draws = selected_draws))
  } else {
    return(list(mat = mat, draws = sub_df$posterior))
  }
}

process_metso_slice_averaged <- function(metso_sub_ids, metso_precomputed, control_precomputed, g_id, n_post, draw_names, mode_comp) {
  metso_results <- list()
  counter <- 1
  
  for (m_id in metso_sub_ids) {
    prep_m <- metso_precomputed[[as.character(m_id)]]
    if (is.null(prep_m)) next
    
    pair_collector <- list()
    p_counter <- 1
    
    for (c_id in names(control_precomputed)) {
      prep_c <- control_precomputed[[c_id]]
      if (is.null(prep_c)) next
      
      res_df <- tryCatch({
        calculate_pairwise_metrics_fast(
          metso_prep = prep_m,
          bau_prep = prep_c,
          num_posteriors = n_post
        )
      }, error = function(e) return(NULL))
      
      if (!is.null(res_df)) {
        pair_collector[[p_counter]] <- res_df
        p_counter <- p_counter + 1
      }
    }
    
    if (length(pair_collector) > 0) {
      combined_dt <- data.table::rbindlist(pair_collector)
      
      if (length(draw_names) == n_post && !identical(draw_names, "mean")) {
        draw_map <- data.table(posterior = seq_len(n_post), actual_draw = draw_names)
        combined_dt <- merge(combined_dt, draw_map, by = "posterior")
        combined_dt[, posterior := actual_draw][, actual_draw := NULL]
      }
      
      # Calculate valid pairwise comparisons
      # When n_valid < 2, sample variance is mathematically undefined -> assign NA_real_
      # When n_valid >= 2, apply a floor (1e-4) to genuine zero variance
      summary_dt <- combined_dt[, {
        valid_vals <- val[!is.na(val)]
        n_valid    <- length(valid_vals)
        
        list(
          val_mean   = if (n_valid >= 1) mean(valid_vals) else NA_real_,
          val_sd     = if (n_valid >= 2) pmax(sd(valid_vals), 1e-4) else NA_real_,
          n_controls = n_valid
        )
      }, by = .(posterior, metrics)]
      
      # Filter records according to mode:
      # - Mode "indomain_both": strictly require valid sd and >= 2 controls
      # - Mode "full": retain rows, preserving NA_real_ for single-control cases
      if (mode_comp == "indomain_both") {
        summary_dt <- summary_dt[!is.na(val_sd) & n_controls >= 2]
      } else {
        summary_dt <- summary_dt[!is.na(val_mean)]
      }
      
      if (nrow(summary_dt) > 0) {
        summary_dt[, standid_treated := m_id]
        summary_dt[, group_number := g_id]
        
        metso_results[[counter]] <- summary_dt
        counter <- counter + 1
      }
    }
  }
  
  if (length(metso_results) > 0) return(data.table::rbindlist(metso_results))
  return(NULL)
}


# 4. Processing Loop
cols_drop <- c("scenario", "standid", "posterior", "year", "reg", "trees")

for (g_id in common_groups) {
  out_file <- file.path(output_dir, paste0("metrics_group_", g_id, ".parquet"))
  if (file.exists(out_file)) {
    cat(sprintf("Group %s already processed. Skipping.\n", g_id))
    next
  }
  
  group_stands <- groups_df[groups_df$group_number == g_id, ]
  metso_ids <- group_stands$standid[group_stands$metso == 1]
  control_ids <- group_stands$standid[group_stands$metso == 0]
  if (length(metso_ids) == 0 || length(control_ids) == 0) next
  
  g_year  <- group_stands$year[1]
  g_reg   <- group_stands$reg[1]
  g_trees <- group_stands$trees[1]
  
  group_preds <- predY_ds |> 
    filter(
      year == g_year,
      reg == g_reg,
      trees == g_trees,
      standid %in% group_stands$standid
    ) |> 
    select(scenario, standid, posterior, all_of(name_spp)) |> 
    collect()
  
  group_pa <- predPA_ds |> 
    filter(
      year == g_year,
      reg == g_reg,
      trees == g_trees,
      standid %in% group_stands$standid
    ) |> 
    select(scenario, standid, posterior, all_of(name_spp)) |> 
    collect()
  
  if (run_mode == "draws") {
    group_preds <- group_preds |> filter(posterior %in% selected_draws)
    group_pa    <- group_pa |> filter(posterior %in% selected_draws)
  }
  
  group_preds_dt <- data.table::as.data.table(group_preds)
  data.table::setkeyv(group_preds_dt, c("scenario", "standid", "posterior"))
  
  group_pa_dt <- data.table::as.data.table(group_pa)
  data.table::setkeyv(group_pa_dt, c("scenario", "standid", "posterior"))
  
  stands_manifest <- rbind(
    data.table::data.table(standid = metso_ids, scenario = "metso"),
    data.table::data.table(standid = control_ids, scenario = "control")
  )
  
  stand_items <- split(stands_manifest, seq_len(nrow(stands_manifest)))
  
  precomputed_results <- parallel::mclapply(
    X = stand_items,
    FUN = function(item) {
      st_id <- item$standid
      scen  <- item$scenario
      
      p_sub  <- group_preds_dt[.(scen, st_id), nomatch = NULL]
      pa_sub <- group_pa_dt[.(scen, st_id), nomatch = NULL]
      if (nrow(p_sub) == 0 || nrow(pa_sub) == 0) return(NULL)
      
      prep_hurdle <- prepare_stand_matrix(as.data.frame(p_sub), run_mode, n_draws, cols_drop, selected_draws)
      prep_pa     <- prepare_stand_matrix(as.data.frame(pa_sub), run_mode, n_draws, cols_drop, selected_draws)
      metrics     <- precalculate_stand_metrics(prep_hurdle$mat, prep_pa$mat, TrData_processed, divmax)
      
      list(
        standid = as.character(st_id),
        scenario = scen,
        metrics = metrics
      )
    },
    mc.cores = n_cores
  )
  
  metso_precomputed <- list()
  control_precomputed <- list()
  
  for (res in precomputed_results) {
    if (is.null(res)) next
    if (res$scenario == "metso") {
      metso_precomputed[[res$standid]] <- res$metrics
    } else {
      control_precomputed[[res$standid]] <- res$metrics
    }
  }
  
  valid_metso_ids <- names(metso_precomputed)
  valid_control_ids <- names(control_precomputed)
  
  # Secondary check on valid controls
  min_controls_req <- if (comparison_mode == "indomain_both") 2 else 1
  if (length(valid_metso_ids) == 0 || length(valid_control_ids) < min_controls_req) {
    cat(sprintf("Group %s lacks sufficient evaluated stands across scenarios. Skipping.\n", g_id))
    next
  }
  
  n_post <- nrow(metso_precomputed[[1]]$pred_mat)
  draw_names <- if (run_mode == "mean") "mean" else if (run_mode == "draws") selected_draws else seq_len(n_post)
  
  chunk_size <- ceiling(length(valid_metso_ids) / n_cores)
  metso_slices <- split(valid_metso_ids, ceiling(seq_along(valid_metso_ids) / chunk_size))
  
  results <- parallel::mclapply(
    X = metso_slices,
    FUN = function(slice) {
      process_metso_slice_averaged(
        metso_sub_ids       = slice,
        metso_precomputed   = metso_precomputed,
        control_precomputed = control_precomputed,
        g_id                = g_id,
        n_post              = n_post,
        draw_names          = draw_names,
        mode_comp           = comparison_mode
      )
    },
    mc.cores = n_cores
  )
  
  results_dt <- data.table::rbindlist(results)
  if (nrow(results_dt) > 0) {
    arrow::write_parquet(results_dt, out_file)
    cat(sprintf("[%s] Group %s written: %d treated stands evaluated against %d controls.\n", 
                format(Sys.time(), "%H:%M:%S"), g_id, length(valid_metso_ids), length(valid_control_ids)))
  }
  
  rm(group_preds, group_pa, group_preds_dt, group_pa_dt, stands_manifest, stand_items, precomputed_results, 
     metso_precomputed, control_precomputed, results, results_dt)
  gc()
}

message("    End time: ", Sys.time())