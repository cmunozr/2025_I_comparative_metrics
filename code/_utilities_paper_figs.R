library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(readxl)
library(scales)
library(dplyr)

# -------------------------------------------------------------------------
# 1. Variance Partitioning
# -------------------------------------------------------------------------

format_vp_csv <- function(file_path, model_type_label) {
  df_raw <- read.csv(file_path, row.names = 1, check.names = FALSE)
  df_filtered <- df_raw[!rownames(df_raw) %in% c("Intercept", "R2", "r2", "AUC", "RMSE", "TjurR2"), ]
  
  df_filtered$Covariate_Group <- rownames(df_filtered)
  
  df_long <- df_filtered |>
    tidyr::pivot_longer(
      cols = -Covariate_Group,
      names_to = "Species",
      values_to = "Variance_Prop"
    ) |>
    dplyr::mutate(
      Model_Component = model_type_label,
      Effect_Type = dplyr::case_when(
        stringr::str_detect(tolower(Covariate_Group), "random|vakio|site|route|year|sample") ~ "Random",
        TRUE ~ "Fixed"
      ),
      Covariate_Clean = dplyr::case_when(
        Covariate_Group == "forest_structure" ~ "Forest Structure",
        Covariate_Group == "climate"          ~ "Climate",
        Covariate_Group == "survey"           ~ "Survey Effort",
        Covariate_Group == "Random: vakio"    ~ "Route (Random)",
        Covariate_Group == "Random: year"     ~ "Year (Random)",
        TRUE ~ Covariate_Group
      )
    )
  return(df_long)
}

create_vp_boxplot_panel <- function(data_subset, 
                                    panel_title, 
                                    fixed_color = "#E69F00", 
                                    random_color = "#0072B2", 
                                    group_order = NULL) {
  
  if (!is.null(group_order)) {
    data_subset$Covariate_Clean <- factor(data_subset$Covariate_Clean, levels = group_order)
  }
  
  summary_means <- data_subset |>
    dplyr::group_by(Covariate_Clean, Effect_Type) |>
    dplyr::summarize(mean_val = mean(Variance_Prop, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(label_y = -0.04, label_text = sprintf("%.2f", mean_val))
  
  fixed_count <- length(unique(data_subset$Covariate_Clean[data_subset$Effect_Type == "Fixed"]))
  divider_pos <- fixed_count + 0.5
  
  p <- ggplot(data_subset, aes(x = Covariate_Clean, y = Variance_Prop)) +
    geom_boxplot(aes(fill = Effect_Type), width = 0.38, outlier.shape = NA, alpha = 0.85, color = "#222222", linewidth = 0.35) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.45, color = "#333333") +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 1.8, fill = "white", color = "black", stroke = 0.6) +
    geom_text(data = summary_means, aes(x = Covariate_Clean, y = label_y, label = label_text), 
              size = 2.6, inherit.aes = FALSE, fontface = "bold", family = "sans") +
    geom_vline(xintercept = divider_pos, linetype = "dashed", color = "#666666", linewidth = 0.4) +
    scale_fill_manual(values = c("Fixed" = fixed_color, "Random" = random_color)) +
    scale_y_continuous(limits = c(-0.06, 1.01), breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = panel_title, 
      x = NULL, 
      y = "Proportion of Explained Variation"
    ) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b = 6)),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "black", size = 7.5),
      axis.text.y = element_text(color = "black", size = 7.5),
      axis.title.y = element_text(size = 8, margin = margin(r = 6)),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 8, l = 10)
    )
  return(p)
}

# File Paths
pa_file <- "models/fbs_M016PA_thin_150_samples_1000_chains_4/parameter_estimates/parameter_estimates_VP.csv"
ab_file <- "models/fbs_M016_thin_250_samples_1000_chains_4/parameter_estimates/parameter_estimates_VP.csv"

df_pa <- format_vp_csv(pa_file, "Presence-Absence")
df_ab <- format_vp_csv(ab_file, "Conditional Abundance")
covariate_order <- c("Forest Structure", "Climate", "Survey Effort", "Route (Random)", "Year (Random)")

# Generate Panels with titles aligned to the plot body
p_a <- create_vp_boxplot_panel(df_pa, "(a) Presence-Absence", group_order = covariate_order)
p_b <- create_vp_boxplot_panel(df_ab, "(b) Conditional Abundance", group_order = covariate_order)

# Combine panels
two_panel_vp <- p_a | p_b

output_dir <- file.path("models", "hurdle_summary_figures")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Export at Elsevier 190 mm width
ggsave(
  filename = file.path(output_dir, "figure_variance_partitioning_hurdle.pdf"),
  plot = two_panel_vp,
  width = 190,
  height = 95,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "figure_variance_partitioning_hurdle.png"),
  plot = two_panel_vp,
  width = 190,
  height = 95,
  units = "mm",
  dpi = 500
)

# -------------------------------------------------------------------------
# 2. Beta figure
# -------------------------------------------------------------------------

pa_file <- "models/fbs_M016PA_thin_150_samples_1000_chains_4/parameter_estimates/parameter_estimates_Beta.xlsx"
pa_cov_keys <- "models/fbs_M016PA_thin_150_samples_1000_chains_4/parameter_estimates/covariate_key_table.csv"

ab_file <- "models/fbs_M016_thin_250_samples_1000_chains_4/parameter_estimates/parameter_estimates_Beta.xlsx"
ab_cov_keys <- "models/fbs_M016_thin_250_samples_1000_chains_4/parameter_estimates/covariate_key_table.csv"

spp_data <- "models/diagnostic_models_df.csv"

# Helper function to abbreviate scientific names and resolve colliding binomials
disambiguate_species <- function(spp_vector) {
  parts <- str_split(spp_vector, "\\s+", simplify = TRUE)
  gen <- parts[, 1]
  spe <- parts[, 2]
  
  # Standard 1-letter abbreviation
  short_standard <- paste0(substr(gen, 1, 1), ". ", spe)
  
  # Identify collision cases (e.g., Poecile montanus vs Passer montanus)
  dup_abbr <- short_standard[duplicated(short_standard)]
  
  # Apply 3-letter genus abbreviation specifically to colliding taxa
  case_when(
    short_standard %in% dup_abbr ~ paste0(substr(gen, 1, 3), ". ", spe),
    TRUE ~ short_standard
  )
}

# Prepare species lookup table from diagnostic data
df_spp_traits <- read.csv(spp_data, stringsAsFactors = FALSE) %>%
  select(species, Specialist) %>%
  distinct() %>%
  mutate(
    species_short = disambiguate_species(species),
    Specialist_group = factor(Specialist, levels = c("Generalist", "Specialist"))
  )

# Ordered factor levels: Generalists on top, Specialists below; A-Z top to bottom
ordered_levels <- df_spp_traits %>%
  arrange(Specialist_group, desc(species_short)) %>%
  pull(species_short)

# Data Preparation Function
prepare_beta_data <- function(beta_path, key_path, part_label, spp_lookup, factor_levels) {
  df_beta <- readxl::read_excel(beta_path)
  df_keys <- read.csv(key_path, stringsAsFactors = FALSE)
  
  cov_col <- intersect(c("Original_Covariate", "Model_Covariate", "model_cov"), colnames(df_keys))[1]
  idx_col <- intersect(c("Plot_Index", "index", "Index"), colnames(df_keys))[1]
  inline_col <- intersect(c("Inline", "inline", "Inline_Label", "label"), colnames(df_keys))[1]
  
  if (is.na(cov_col) || is.na(idx_col)) {
    stop(paste("Required columns not identified in:", key_path))
  }
  
  # If Inline column is present, use it; otherwise fallback to Plot_Index
  if (!is.na(inline_col)) {
    df_keys <- df_keys %>%
      rename(Original_Covariate = all_of(cov_col),
             Plot_Index = all_of(idx_col),
             Inline_Label = all_of(inline_col))
  } else {
    df_keys <- df_keys %>%
      rename(Original_Covariate = all_of(cov_col),
             Plot_Index = all_of(idx_col)) %>%
      mutate(Inline_Label = Plot_Index)
  }
  
  # Define conceptual covariate categories and ordering
  cov_group_levels <- c("G1", "G2", "G3", "G4")
  
  df_keys <- df_keys %>%
    filter(Original_Covariate != "(Intercept)") %>%
    mutate(
      Covariate_Group = case_when(
        str_detect(Original_Covariate, "diameter|stand_age") ~ "G1",
        str_detect(Original_Covariate, "broadleaves|volume_spruce") ~ "G2",
        str_detect(Original_Covariate, "canopy_cover_whole|height|tree_extent|stand_length") ~ "G3",
        str_detect(Original_Covariate, "tree_remove") ~ "G4",
        TRUE ~ "Other"
      ),
      Covariate_Group = factor(Covariate_Group, levels = cov_group_levels)
    )
  
  # Determine ordered factor levels for the x-axis
  ordered_cov_labels <- df_keys %>%
    arrange(Covariate_Group, as.numeric(str_extract(Plot_Index, "\\d+"))) %>%
    pull(Inline_Label) %>%
    unique()
  
  df_beta %>%
    inner_join(df_keys, by = c("model_cov" = "Original_Covariate")) %>%
    left_join(spp_lookup, by = c("Species" = "species")) %>%
    filter(type != "Intercept") %>%
    mutate(
      model_part = part_label,
      support_category = case_when(
        support >= 0.80 ~ "+ Positive",
        supportNeg >= 0.80 ~ "- Negative",
        TRUE ~ "No Support"
      ),
      support_category = factor(
        support_category,
        levels = c("+ Positive", "No Support", "- Negative")
      ),
      Inline_Label = factor(Inline_Label, levels = ordered_cov_labels),
      species_short = factor(species_short, levels = factor_levels)
    )
}

# Load and process data
df_pa <- prepare_beta_data(pa_file, pa_cov_keys, "Presence-Absence", df_spp_traits, ordered_levels)
df_ab <- prepare_beta_data(ab_file, ab_cov_keys, "Conditional Abundance", df_spp_traits, ordered_levels)

# High-contrast tonal palette matching the sign-support tiers
support_palette <- c(
  "+ Positive" = "#E69F00",
  "No Support" = "#e5e5e5",
  "- Negative" = "#0072B2"
)

# Custom Elsevier journal theme
theme_elsevier_beta <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.01),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.01),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(face = "italic", size = 7, color = "black"),
    axis.title = element_text(size = 9, face = "plain", color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(b = 5)),
    strip.placement = "outside",
    strip.background = element_rect(fill = "gray96", color = "black", linewidth = 0.04),
    strip.text.y = element_text(size = 8.5, face = "bold", angle = 90),
    strip.text.x = element_text(size = 8, face = "bold"),
    panel.spacing.x = unit(1.5, "pt"),
    panel.spacing.y = unit(2, "pt"),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3.5, "mm"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 2)
  )

# Construct Biplot Panels
p_pa <- ggplot(df_pa, aes(x = Inline_Label, y = species_short, fill = support_category)) +
  geom_tile(color = "gray92", linewidth = 0.15, width = 1, height = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(Specialist_group ~ Covariate_Group, scales = "free", space = "free", switch = "y") +
  scale_fill_manual(values = support_palette, drop = FALSE, name = "Forest-structure\neffect") +
  labs(title = "(a) Presence-Absence", x = NULL, y = NULL) +
  theme_elsevier_beta +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_ab <- ggplot(df_ab, aes(x = Inline_Label, y = species_short, fill = support_category)) +
  geom_tile(color = "gray92", linewidth = 0.15, width = 1, height = 1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(Specialist_group ~ Covariate_Group, scales = "free", space = "free") +
  scale_fill_manual(values = support_palette, drop = FALSE, name = "Posterior\nSupport") +
  labs(title = "(b) Conditional Abundance", x = NULL, y = NULL) +
  theme_elsevier_beta +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

# Combine panels with shared layout and legend
biplot_beta_hurdle <- (p_pa + p_ab) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Wrapper with unified bottom x-axis title
final_figure <- wrap_elements(biplot_beta_hurdle) +
  labs(tag = "Forest-structure covariates") +
  theme(
    plot.tag.position = "bottom",
    plot.tag = element_text(size = 9, face = "plain", hjust = 0.45, margin = margin(t = 2, b = 2))
  )

# Export
ggsave(
  filename = "models/hurdle_summary_figures/figure_beta_estimates_hurdle.pdf",
  plot = final_figure,
  width = 190,
  height = 230,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = "models/hurdle_summary_figures/figure_beta_estimates_hurdle.png",
  plot = final_figure,
  width = 190,
  height = 230,
  units = "mm",
  dpi = 500
)

# -------------------------------------------------------------------------
# 3. consistency heatmaps
# -------------------------------------------------------------------------

footprints <- c("E_MSA_loss", "E_PDF_unforced", "E_PDF_forced")
ecological <- c("EG_loss", "E_FRic_loss", "E_RaoQ_loss", "E_BCD")
metrics    <- c(footprints, ecological)

overall_rank_summary <- read.csv(file.path("results","consistency","q2_rank_summary.csv"))
stand_summary <- read.csv(file.path("results","consistency","q2_stand_discordance_summary.csv"))

output_dir <- file.path("results", "consistency")


# --- 5. Visualization Heatmaps (Elsevier / Ecological Indicators Style) ---
cat("\n--- Generating Heatmap Figures ---\n")

# Publication-grade label dictionaries matching the mock-up
footprint_labels <- c(
  "E_MSA_loss"     = "EMSA",
  "E_PDF_unforced" = "EPDF-U",
  "E_PDF_forced"   = "EPDF-F"
)

ecological_labels <- c(
  "EG_loss"     = "EG",
  "E_FRic_loss" = "EFRic",
  "E_RaoQ_loss" = "ERaoQ",
  "E_BCD"       = "EBCD"
)

# Standardized theme matching S06b / Elsevier specifications
theme_elsevier_heatmap <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.text.x        = element_text(size = 8, face = "plain", color = "black"),
    axis.text.y        = element_text(size = 8, face = "plain", color = "black"),
    axis.title         = element_blank(),
    axis.ticks         = element_blank(),
    plot.title         = element_text(size = 8.5, face = "bold", hjust = 0, margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 7, color = "gray25", hjust = 0, margin = margin(b = 4)),
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.key.width   = unit(3, "mm"),
    legend.key.height  = unit(7, "mm"),
    plot.margin        = margin(t = 6, r = 4, b = 6, l = 4)
  )

# Base datasets filtered and formatted with ordered factor levels
df_rank_fig <- overall_rank_summary |>
  filter(metric_a %in% footprints, metric_b %in% ecological) |>
  mutate(
    metric_a = factor(metric_a, levels = footprints, labels = footprint_labels),
    metric_b = factor(metric_b, levels = rev(ecological), labels = rev(ecological_labels))
  )

df_sign_fig <- stand_summary |>
  filter(metric_a %in% footprints, metric_b %in% ecological) |>
  mutate(
    metric_a = factor(metric_a, levels = footprints, labels = footprint_labels),
    metric_b = factor(metric_b, levels = rev(ecological), labels = rev(ecological_labels))
  ) |>
  tidyr::complete(
    metric_a = factor(footprint_labels, levels = footprint_labels),
    metric_b = factor(rev(ecological_labels), levels = rev(ecological_labels))
  )

# --- Panel (a): Median Spearman Rank Correlation ---
p_a <- ggplot(df_rank_fig, aes(x = metric_a, y = metric_b, fill = median_rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", median_rho)), size = 2.7, family = "sans") +
  scale_fill_gradient2(
    limits = c(-1, 1), midpoint = 0,
    low = "#b2182b", mid = "#eef0f6", high = "#386cb0",
    na.value = "grey90"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), drop = FALSE) +
  labs(
    title = "(a) Rank concordance",
    subtitle = expression(Median~Spearman~rho),
    fill = expression(rho)
  ) +
  theme_elsevier_heatmap

# --- Panel (b): Median Top-10% Overlap Lift ---
p_b <- ggplot(df_rank_fig, aes(x = metric_a, y = metric_b, fill = median_top_lift)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", median_top_lift)), size = 2.7, family = "sans") +
  scale_fill_gradient2(
    midpoint = 1.0,
    low = "#b2182b", mid = "#f5f5f5", high = "#438473",
    na.value = "grey90"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), drop = FALSE) +
  labs(
    title = "(b) Top-10% overlap",
    subtitle = "Overlap lift vs. chance",
    fill = "Lift"
  ) +
  theme_elsevier_heatmap +
  theme(axis.text.y = element_blank())

# --- Composite Figure
fig_consistency_3panel <- p_a + p_b +
  plot_layout(ncol = 2, widths = c(1.15, 1, 1))

ggsave(
  filename = file.path(output_dir, "figure_q2_metric_consistency_3panel.pdf"),
  plot     = fig_consistency_3panel,
  width    = 190,
  height   = 75,
  units    = "mm",
  device   = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "figure_q2_metric_consistency_3panel.png"),
  plot     = fig_consistency_3panel,
  width    = 190,
  height   = 75,
  units    = "mm",
  dpi      = 500
)
