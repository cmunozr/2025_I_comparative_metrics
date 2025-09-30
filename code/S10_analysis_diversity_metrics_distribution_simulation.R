# Install packages if you don't have them
# install.packages("tidyverse")
# install.packages("viridis")

library(tidyverse) # Includes ggplot2, dplyr, tidyr, etc.
library(viridis)   # For the color palette

# --- Data Simulation ---
n_posterior_samples <- 1000  # Number of posterior samples
set.seed(11072024) # For reproducibility

metrics <- c('MSA', 'Sorensen', 'Geometric', 'Rao', 'BetaDiv')
sampling_conditions <- c('Complete', 'Partial') # 'Complete' likely means full/reference sampling

# We will store the partial data frames here
results_list <- list()

# Define the characteristics of the simulated distributions
# (mean, standard deviation) for each metric and condition
# Using a nested list in R
simulation_params <- list(
  'MSA' = list('Complete' = c(mean = 0.8, sd = 0.2), 'Partial' = c(mean = 0.3, sd = 0.3)),
  'Sorensen'  = list('Complete' = c(mean = 1.0, sd = 0.15), 'Partial' = c(mean = 0.9, sd = 0.2)),
  'Geometric'  = list('Complete' = c(mean = 0.4, sd = 0.25), 'Partial' = c(mean = 0.1, sd = 0.35)),
  'Rao'  = list('Complete' = c(mean = 0.1, sd = 0.3), 'Partial' = c(mean = 0.05, sd = 0.35)),
  'BetaDiv' = list('Complete' = c(mean = 0.6, sd = 0.18), 'Partial' = c(mean = 0.58, sd = 0.19))
)

# Iterate to generate simulated data
for (metric in metrics) {
  for (condition in sampling_conditions) {
    params <- simulation_params[[metric]][[condition]]
    mean_val <- params['mean']
    sd_val <- params['sd']
    
    # Simulate the posterior distribution of the average effect size
    effect_samples <- rnorm(n_posterior_samples, mean = mean_val, sd = sd_val)
    
    # Create a temporary data frame and add it to the list
    temp_df <- data.frame(
      EffectSize = effect_samples,
      Metric = metric,
      Sampling = condition
    )
    results_list[[length(results_list) + 1]] <- temp_df
  }
}

# Combine all results into a single DataFrame
all_results_df <- bind_rows(results_list) %>%
  # Convert to factors to control order and labels in ggplot
  mutate(
    Metric = factor(Metric, levels = metrics),
    Sampling = factor(Sampling, levels = sampling_conditions)
  )


# --- Creating the Plot with ggplot2 ---

p <- ggplot(all_results_df, aes(x = Metric, y = EffectSize, fill = Sampling)) +
  # Use geom_violin
  # position_dodge places violins side-by-side
  # scale="width" makes all violins have the same maximum width
  # draw_quantiles adds lines for quartiles (0.25, 0.50, 0.75)
  geom_violin(position = position_dodge(width = 0.9), scale = "width", trim = FALSE,
              draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.8) +
  
  # Add a horizontal line at y=0 for reference (no effect)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  
  # Use the Viridis color palette (similar to the Python version)
  scale_fill_(name = "Sampling\nIntensity") + # Name for the legend
  
  # Titles and labels
  labs(
    title = "Comparison of Mean Effect Size (METSO - Bussiness as usual)",
    subtitle = "Under Different Sampling Intensities (Local Scale)",
    x = "Biodiversity Metric",
    y = "Average Effect Size\n(Scenario A - Baseline)"
  ) +
  
  # A clean theme for presentation
  theme_minimal(base_size = 14) +
  # Rotate x-axis labels if necessary
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "C:/Users/Carlos Munoz/Documents/Ph.D/6_courses/_meetings/2025/five_figure.jpeg",    # The name of the file you want to save
  plot = p,           # The plot object you want to save
  width = 7,          # The width of the plot in inches
  height = 5,         # The height of the plot in inches
  units = "in",       # The units for the width and height
  dpi = 300           # The resolution in dots per inch (PPI)
)

# --- Interpretation (Example for your presentation) ---
# When presenting this graph, you could say something like:
# - "Here we see the simulated posterior distributions of the average effect of Scenario A for 5 metrics, comparing full ('Complete') sampling with reduced ('Partial') sampling."
# - "The dashed line at zero indicates no effect. We look for distributions clearly above or below zero."
# - "We can observe that metrics like 'Shannon' and 'PhyloDiv' appear robust, detecting a clear positive effect even with partial sampling (green and blue violins are similar and far from zero)."
# - "On the other hand, 'Richness' and 'Simpson' show a strong signal with complete data (Complete), but their detection ability degrades noticeably (the distribution moves closer to zero and widens) with partial sampling."
# - "'FuncDiv', in this simulation, seems insensitive to the change induced by Scenario A, as its distributions are centered near zero in both cases."
# - "This type of plot allows us to visually assess which metrics are more sensitive and robust to changes in sampling effort for detecting the impact of a specific scenario."