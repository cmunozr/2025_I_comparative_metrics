# Install packages if you don't have them
# install.packages("tidyverse")
# install.packages("viridis")

library(tidyverse) # Includes ggplot2, dplyr, tidyr, etc.
library(viridis)   # For the color palette
library(ggplot2)
library(here)

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

set.seed(11072024)
test <- T

if(test){
  num_loc <- 100
  samples <- 100
  to_sample <- list("a" = sample(seq(1971), size = num_loc), "b"= sample(seq(4000), size = samples))
}


# --- Data ---

metric_files <- list.files(file.path("results", "metrics"), pattern = ".rds", full.names = T)
metric_files <- metric_files[stringr::str_detect(metric_files, pattern = "TrData", negate = T)]
n_sample <- 1971

metrics <- basename(metric_files) |> 
  gsub(pattern = ".rds", replacement = "") |> 
  toupper()

sampling_conditions <- c("local") # 'Complete' means full/reference sampling

# We will store the partial data frames here
results_list <- lapply(X = metric_files, 
                       FUN = function(X, is.test = test, to_s = to_sample){
                           effect <- readRDS(X)[["norm"]] 
                           if(is.test) effect <- effect[to_s$a, to_s$b] 
                           if(is.data.frame(effect)) effect <- as.matrix(effect)
                           effect <- as.vector(effect)
                           if(length(effect) != 10000){
                             return(NA)# hardcode for example
                           }else{
                             res <- data.frame(EffectSize = effect)
                             return(res)
                             } 
                           }
                       )
  
# Combine all results into a single DataFrame
all_results_df <- bind_rows(results_list) %>%
  # Convert to factors to control order and labels in ggplot
  mutate(
    Metric = factor(rep(x = metrics, each = 10000), levels = metrics),
    Sampling = rep(x = sampling_conditions, 60000)
  )

all_results_df <- all_results_df[!is.infinite(all_results_df$EffectSize), ]

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
  scale_fill_discrete(name = "Sampling\nIntensity") + # Name for the legend
  
  # Titles and labels
  labs(
    title = "Comparison of Mean Effect Size (Bussiness as usual - METSO)",
    subtitle = "Under Different Sampling Intensities (Local Scale)",
    x = "Biodiversity Metric",
    y = "Average Effect Size\n(BAU - Metso)"
  ) +
  
  # A clean theme for presentation
  theme_minimal(base_size = 14) +
  # Rotate x-axis labels if necessary
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(file.path("analysis"), showWarnings = F)

ggsave(
  "analysis/comparisson_base_figure.jpeg",    # The name of the file you want to save
  plot = p,           # The plot object you want to save
  width = 7,          # The width of the plot in inches
  height = 5,         # The height of the plot in inches
  units = "in",       # The units for the width and height
  dpi = 300           # The resolution in dots per inch (PPI)
)

# --- Interpretation ---
# - "Here we see the simulated posterior distributions of the average effect of Scenario A for 5 metrics, comparing full ('Complete') sampling with reduced ('Partial') sampling."
# - "The dashed line at zero indicates no effect. We look for distributions clearly above or below zero."
# - "We can observe that metrics like 'Shannon' and 'PhyloDiv' appear robust, detecting a clear positive effect even with partial sampling (green and blue violins are similar and far from zero)."
# - "On the other hand, 'Richness' and 'Simpson' show a strong signal with complete data (Complete), but their detection ability degrades noticeably (the distribution moves closer to zero and widens) with partial sampling."
# - "'FuncDiv', in this simulation, seems insensitive to the change induced by Scenario A, as its distributions are centered near zero in both cases."
# - "This type of plot allows us to visually assess which metrics are more sensitive and robust to changes in sampling effort for detecting the impact of a specific scenario."