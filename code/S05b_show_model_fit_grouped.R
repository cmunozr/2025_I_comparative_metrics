# --- 1. Load Libraries ---
library(colorspace)
library(tidyverse) 
library(gridExtra) 
library(Hmsc)

# --- 2. Configuration ---
base_dir <- "models"
ensure_dir <- file.path(base_dir)
dir.create(ensure_dir, recursive = TRUE, showWarnings = FALSE)

# Extract species names directly from a reference fitted Hmsc object
reference_fitted_model <- file.path(
  base_dir, 
  "fbs_M016PA_thin_150_samples_1000_chains_4", 
  "fitted_fbs_M016PA_thin_150_samples_1000_chains_4.rds"
)

if (file.exists(reference_fitted_model)) {
  hm_ref <- readRDS(reference_fitted_model)
  species_names <- colnames(hm_ref$Y)
  message("Extracted ", length(species_names), " species names from reference model.")
} else {
  # Fallback: search for any fitted model RDS file in base_dir
  fitted_files <- list.files(base_dir, pattern = "^fitted_.*\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(fitted_files) > 0) {
    hm_ref <- readRDS(fitted_files[1])
    species_names <- colnames(hm_ref$Y)
    message("Extracted ", length(species_names), " species names from: ", basename(fitted_files[1]))
  } else {
    species_names <- NULL
    warning("No fitted Hmsc object found. Defaulting to numeric indices.")
  }
}

# Define cross-validation / hold-out strategies to inspect
validation_strategies <- c("route_blocked_cv", "metso_holdout", "north_south") 

strategy_map <- c(
  "metso_holdout"    = "ho_metso",
  "route_blocked_cv" = "cv_route",
  "north_south"      = "north_south"
)

# --- 3. Folder Processing Function ---
process_model_folder_by_strategy <- function(folder_path, strategies, strat_map, spp_vec = NULL) {
  run_name <- basename(folder_path)
  fit_dir  <- file.path(folder_path, "model_fit")
  
  if (!dir.exists(fit_dir)) return(NULL)
  
  files     <- list.files(fit_dir, full.names = TRUE)
  path_waic <- files[str_detect(files, paste0("waic_", run_name, "\\.rds"))][1]
  path_mf   <- files[str_detect(files, paste0("mf_", run_name, "\\.rds"))][1]
  
  if (is.na(path_waic) || is.na(path_mf)) return(NULL)
  
  short_model_name <- str_extract(run_name, "fbs_M[0-9]{3}(\\.[1-9])?")
  meta_thin        <- str_extract(run_name, "thin_[0-9]+")
  meta_samples     <- str_extract(run_name, "samples_[0-9]+")
  meta_chains      <- str_extract(run_name, "chains_[0-9]+")
  model_class      <- paste(meta_thin, meta_samples, meta_chains, sep = "_")
  
  model_response_type <- if (str_detect(run_name, "PA")) "Presence-Absence" else "Continuous Abundance"
  
  waic_val <- read_rds(path_waic) %>% unlist() %>% as.numeric()
  mf_raw   <- tryCatch(read_rds(path_mf), error = function(e) NULL)
  
  # Determine species vector from local fitted model if available
  local_fitted <- list.files(folder_path, pattern = "^fitted_.*\\.rds$", full.names = TRUE)[1]
  if (!is.na(local_fitted) && file.exists(local_fitted)) {
    hm_local <- readRDS(local_fitted)
    spp_local <- colnames(hm_local$Y)
  } else {
    spp_local <- spp_vec
  }
  
  df_mf <- NULL
  if (!is.null(mf_raw)) {
    n_spp <- length(mf_raw[[1]])
    spp_col <- if (!is.null(spp_local) && length(spp_local) == n_spp) spp_local else paste0("sp_", seq_len(n_spp))
    
    df_mf <- do.call("cbind", mf_raw) %>% 
      as.data.frame() %>%
      mutate(
        Model        = short_model_name,
        Model_Type   = model_response_type,
        WAIC         = waic_val,
        Class        = model_class,
        Type         = "Explanatory",
        Strategy     = "Full_Model",
        Species_Idx  = row_number(),
        Species      = spp_col
      )
  }
  
  df_eval_list <- list()
  for (strat in strategies) {
    lbl <- strat_map[strat]
    path_mfeval <- files[str_detect(files, paste0("mfeval_", run_name, "_", lbl, "_\\.rds"))][1]
    
    if (!is.na(path_mfeval) && file.exists(path_mfeval)) {
      mfeval_raw <- tryCatch(read_rds(path_mfeval), error = function(e) NULL)
      
      if (!is.null(mfeval_raw)) {
        n_spp <- length(mfeval_raw[[1]])
        spp_col <- if (!is.null(spp_local) && length(spp_local) == n_spp) spp_local else paste0("sp_", seq_len(n_spp))
        
        df_strat_eval <- do.call("cbind", mfeval_raw) %>% 
          as.data.frame() %>%
          mutate(
            Model        = short_model_name,
            Model_Type   = model_response_type,
            WAIC         = waic_val,
            Class        = model_class,
            Type         = "Predictive",
            Strategy     = strat,
            Species_Idx  = row_number(),
            Species      = spp_col
          )
        df_eval_list[[strat]] <- df_strat_eval
      }
    }
  }
  
  return(list(mf = df_mf, mfeval = bind_rows(df_eval_list)))
}

# --- 4. Data Aggregation & CSV Export ---
message("--- Starting Strategy-Aware Data Aggregation ---")

model_dirs  <- list.dirs(base_dir, recursive = FALSE) %>% str_subset("fbs_M")
all_results <- purrr::map(model_dirs, ~process_model_folder_by_strategy(.x, validation_strategies, strategy_map, species_names)) %>% compact()

mf_df     <- purrr::map(all_results, "mf") %>% compact() %>% bind_rows()
mfeval_df <- purrr::map(all_results, "mfeval") %>% compact() %>% bind_rows()

master_df <- bind_rows(mf_df, mfeval_df)

if (nrow(master_df) > 0) {
  master_df <- master_df %>% 
    group_by(Model_Type) %>% 
    mutate(WAIC_Label = fct_reorder(as.factor(round(WAIC, 1)), WAIC)) %>% 
    ungroup() %>%
    mutate(Strategy = factor(Strategy, levels = c("Full_Model", "random_cv", "route_blocked_cv", "metso_holdout", "north_south")))
  
  # --- CSV Outputs ---
  message("--- Writing CSV Outputs with Verified Species Names ---")
  pa_master_data <- master_df %>% filter(Model_Type == "Presence-Absence")
  ca_master_data <- master_df %>% filter(Model_Type == "Continuous Abundance")
  
  write_csv(pa_master_data, file.path(base_dir, "model_fit_PA_comparison.csv"))
  write_csv(ca_master_data, file.path(base_dir, "model_fit_AbuCP_comparison.csv"))
  write_csv(master_df, file.path(base_dir, "model_fit_master_all_metrics.csv"))
  
  summary_metrics_table <- master_df %>%
    group_by(Model, Model_Type, Strategy, Type) %>%
    summarise(
      WAIC = mean(WAIC, na.rm = TRUE),
      across(
        any_of(c("AUC", "TjurR2", "SR2", "RMSE")),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE),
          q25    = ~quantile(.x, 0.25, na.rm = TRUE),
          q75    = ~quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_species = n(),
      .groups = "drop"
    )
  write_csv(summary_metrics_table, file.path(base_dir, "model_fit_summary_table.csv"))
}

# --- 5. Panel Builder Object ---
build_metric_panel <- function(data_subset, metric_column, strategy_name, limit_y) {
  if (!metric_column %in% colnames(data_subset)) return(NULL)
  
  strat_data <- data_subset %>% filter(Strategy == strategy_name)
  if (nrow(strat_data) == 0 || all(is.na(strat_data[[metric_column]]))) {
    p_empty <- ggplot() + 
      annotate("text", x = 1, y = 1, label = paste0(strategy_name, "\n[File Not Found / Skipped]"), 
               size = 3.5, color = "gray50", fontface = "italic") + 
      theme_void() +
      labs(title = paste0(metric_column, ": ", strategy_name))
    return(p_empty)
  }
  
  display_title <- case_when(
    strategy_name == "Full_Model"       ~ "Explanatory (Full Fit)",
    strategy_name == "random_cv"        ~ "Predictive: Random CV",
    strategy_name == "route_blocked_cv" ~ "Predictive: Route Blocked",
    strategy_name == "metso_holdout"    ~ "Predictive: METSO Holdout",
    strategy_name == "north_south"      ~ "Predictive: North Holdout",
    TRUE                                ~ strategy_name
  )
  
  p <- ggplot(strat_data, aes(x = WAIC_Label, y = .data[[metric_column]], fill = Model)) +
    geom_boxplot(width = 0.4, alpha = 0.6, outlier.size = 0.6) +
    stat_summary(
      fun = mean, geom = "text", 
      aes(label = paste0("\u03BC=", round(after_stat(y), 3))),
      vjust = -1.0, size = 2.0, color = "black", fontface = "bold"
    ) +
    geom_hline(yintercept = ifelse(metric_column == "AUC", 0.5, 0), linetype = "dashed", color = "gray40") +
    coord_cartesian(ylim = limit_y, clip = "off") +
    labs(
      title = display_title,
      x = "WAIC Order",
      y = metric_column
    ) +
    theme_minimal(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  return(p)
}

get_shared_legend <- function(data_subset) {
  p_leg <- ggplot(data_subset, aes(x = WAIC_Label, y = Species_Idx, fill = Model)) + 
    geom_boxplot() + 
    theme_minimal() + 
    theme(legend.position = "bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 8))
  tmp <- ggplot_gtable(ggplot_build(p_leg))
  leg_idx <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  return(tmp$grobs[[leg_idx]])
}

# =========================================================================
# 6. VISUALIZATIONS: PRESENCE-ABSENCE MODELS
# =========================================================================
pa_master <- master_df %>% filter(Model_Type == "Presence-Absence")

if (nrow(pa_master) > 0) {
  message("Generating metric-focused layout for Presence-Absence...")
  pdf_pa <- file.path(base_dir, "model_fit_PA_comparison.pdf")
  pa_pages <- list()
  
  all_strategies <- c("Full_Model", "route_blocked_cv", "metso_holdout", "north_south")
  shared_legend  <- get_shared_legend(pa_master)
  
  # Page 1: Tjur R2
  p_tjur <- purrr::map(all_strategies, ~build_metric_panel(pa_master, "TjurR2", .x, c(-0.05, 1)))
  pa_pages[[1]] <- grid.arrange(grobs = p_tjur, ncol = 4, bottom = shared_legend,
                                top = "Presence-Absence Models: Cross-Strategy Comparison of Tjur R2")
  
  # Page 2: AUC
  p_auc <- purrr::map(all_strategies, ~build_metric_panel(pa_master, "AUC", .x, c(0.45, 1)))
  pa_pages[[2]] <- grid.arrange(grobs = p_auc, ncol = 4, bottom = shared_legend,
                                top = "Presence-Absence Models: Cross-Strategy Comparison of AUC (Discrimination)")
  
  # Page 3: RMSE
  max_rmse_pa <- max(pa_master$RMSE, na.rm = TRUE)
  limit_y_rmse_pa <- c(0, max_rmse_pa * 1.1)
  
  p_rmse <- purrr::map(all_strategies, ~build_metric_panel(pa_master, "RMSE", .x, limit_y_rmse_pa))
  pa_pages[[3]] <- grid.arrange(grobs = p_rmse, ncol = 4, bottom = shared_legend,
                                top = "Presence-Absence Models: Cross-Strategy Comparison of Absolute RMSE")
  
  marrangeGrob(grobs = pa_pages, ncol = 1, nrow = 1, top = NULL) %>%
    ggsave(filename = pdf_pa, width = 14, height = 5.5)
}

# =========================================================================
# 7. VISUALIZATIONS: CONTINUOUS ABUNDANCE MODELS
# =========================================================================
ca_master <- master_df %>% filter(Model_Type == "Continuous Abundance")

if (nrow(ca_master) > 0) {
  message("Generating metric-focused layout for Continuous Abundance...")
  pdf_ca <- file.path(base_dir, "model_fit_AbuCP_comparison.pdf")
  ca_pages <- list()
  
  all_strategies <- c("Full_Model", "route_blocked_cv", "metso_holdout", "north_south")
  shared_legend  <- get_shared_legend(ca_master)
  
  # Page 1: SR2
  p_sr2 <- purrr::map(all_strategies, ~build_metric_panel(ca_master, "SR2", .x, c(-0.2, 1)))
  ca_pages[[1]] <- grid.arrange(grobs = p_sr2, ncol = 4, bottom = shared_legend,
                                top = "Continuous Abundance Models: Cross-Strategy Comparison of Spearman's R2 (SR2)")
  
  # Page 2: RMSE
  max_rmse_ca <- max(ca_master$RMSE, na.rm = TRUE)
  limit_y_rmse_ca <- c(0, max_rmse_ca * 1.1)
  
  p_rmse_ca <- purrr::map(all_strategies, ~build_metric_panel(ca_master, "RMSE", .x, limit_y_rmse_ca))
  pa_pages_ca <- list()
  ca_pages[[2]] <- grid.arrange(grobs = p_rmse_ca, ncol = 4, bottom = shared_legend,
                                top = "Continuous Abundance Models: Cross-Strategy Comparison of Absolute RMSE")
  
  marrangeGrob(grobs = ca_pages, ncol = 1, nrow = 1, top = NULL) %>%
    ggsave(filename = pdf_ca, width = 14, height = 5.5)
}