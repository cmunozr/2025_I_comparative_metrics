library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

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
                                    fixed_color = "#E06666", 
                                    random_color = "#2AA198", 
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