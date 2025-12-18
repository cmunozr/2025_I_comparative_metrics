library(tidyverse)
library(ggtext) 

# 1. Extract and Combine Data
# loop through 'res' list and extract the summary tibble

summary_combined <- map_dfr(names(res), function(nm) {
  
  item <- res[[nm]]
  df <- item[["overall"]]
  df |> 
    as_tibble() |>
    mutate(metric = nm) |>
    select(metric, mean, lower_ci, upper_ci)
})

# 2. Add Significance Flag
# Visually identify if the metric detects a 'non-zero' effect
summary_combined <- summary_combined |>
  mutate(
    significant = ifelse(lower_ci > 0 | upper_ci < 0, "Significant", "Not Significant"),
    # Create a label for the plot
    label = sprintf("%.2f [%.2f; %.2f]", mean, lower_ci, upper_ci)
  )

# 3. Create the Comparative Graphic

p_compare <- ggplot(summary_combined, aes(x = mean, y = reorder(metric, mean))) +
  
  # A. Reference Line (No effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = significant), 
                 height = 0.3, size = 0.8) +
  geom_point(aes(color = significant), size = 4, shape = 18) +
    geom_text(aes(label = label), vjust = -1.5, size = 3, color = "gray30") +
  
  scale_color_manual(values = c("Not Significant" = "gray60", "Significant" = "navyblue")) +
  
  # E. Aesthetics
  labs(
    title = "Comparative Sensitivity of Biodiversity Metrics (Local scale)",
    subtitle = "Posterior Mean Effect Size (METSO - BAU) with 95% CrI",
    x = "Effect Size",
    y = NULL,
    caption = "Note: 'Significant' implies 95% CrI excludes zero."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11)
  )

# 4. Display and Save
print(p_compare)

ggsave(file.path("analysis", "Comparative_Metrics_Summary.png"), 
       plot = p_compare, width = 10, height = 8, bg = "white")