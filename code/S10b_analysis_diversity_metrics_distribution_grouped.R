
library(tidyverse) 
library(viridis)
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

sampling_conditions <- c("local")

# store partial data frames here
results_list <- lapply(X = metric_files, 
                       FUN = function(X, is.test = test, to_s = to_sample){
                           effect <- readRDS(X)[["delta"]] 
                           if(is.test) effect <- effect[to_s$a, to_s$b] 
                           if(is.data.frame(effect)) effect <- as.matrix(effect)
                           effect <- as.vector(effect)
                           if(length(effect) != 10000){
                             return(NA)
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
  # position_dodge places violins side-by-side
  # scale="width" makes all violins have the same maximum width
  # draw_quantiles adds lines for quartiles (0.25, 0.50, 0.75)
  geom_violin(position = position_dodge(width = 0.9), scale = "width", trim = FALSE,
              draw_quantiles = c(0.25, 0.5, 0.75), alpha=0.8) +
  
  # Add a horizontal line at y=0 for reference (no effect)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  
  # Use the Viridis color palette
  scale_fill_discrete(name = "Sampling\nIntensity") + # Name for the legend
  
  # Titles and labels
  labs(
    title = "Comparison of Mean Effect Size (Bussiness as usual - METSO)",
    subtitle = "Under Different Sampling Intensities (Local Scale)",
    x = "Biodiversity Metric",
    y = "Average Effect Size\n(BAU - Metso)"
  ) +
  theme_minimal(base_size = 14) +
  # Rotate x-axis labels if necessary
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(file.path("analysis"), showWarnings = F)

ggsave(
  "analysis/comparisson_base_figure.jpeg",    
  plot = p,           
  width = 7,          
  height = 5,         
  units = "in",       
  dpi = 300          
)
