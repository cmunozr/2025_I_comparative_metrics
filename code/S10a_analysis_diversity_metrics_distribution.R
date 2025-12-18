library(matrixStats) # For efficient row/column summaries
library(tidyverse)
library(glue)
library(gt)

set.seed(11072024)

# --- 1. Call Effect Size (Difference) Posterior Matrix ---

metric_files <- list.files(file.path("results", "metrics"), pattern = ".rds", full.names = T)
metric_files <- metric_files[stringr::str_detect(metric_files, pattern = "TrData", negate = T)]
metric_nms <- basename(metric_files) |> 
  gsub(pattern = ".rds", replacement = "") |> 
  toupper()

metric_nms[which(metric_nms == "PDF")] <- "L-PDF*" 

res <- list()

for(i in 1:length(metric_files)){
  
  metric_nm <- metric_nms[i]
  
  metric_list <- readRDS(metric_files[i])
  
  # This matrix holds the posterior distribution of the effect for each locality
  # here we call the normalized difference
  
  if((metric_nm %in% c("MSA", "SORENSEN","L-PDF*"))){
    effect_size_posterior <- metric_list[["norm"]] |> as.matrix()  
    if(metric_nm == "L-PDF*"){
       effect_size_posterior <- (effect_size_posterior*-1)
    }
  }else{
    effect_size_posterior <- metric_list[["norm"]] |> as.matrix()
  }
  
  effect_size_posterior[is.infinite(effect_size_posterior)] <- NA 
  n_localities <- nrow(effect_size_posterior)
  n_samples <- ncol(effect_size_posterior)
  
  # --- 2. Locality-Level Summaries (for Forest Plot visualization) ---
  
  locality_mean <- rowMeans(effect_size_posterior, na.rm = T)
  locality_ci_lower <- rowQuantiles(effect_size_posterior, probs = 0.025, na.rm = T)
  locality_ci_upper <- rowQuantiles(effect_size_posterior, probs = 0.975, na.rm = T)
  
  # --- 3. Calculate Posterior Distribution of the Overall Average Effect ---
  
  # Average across localities *for each posterior sample*
  # This directly gives the posterior distribution of the mean effect across all localities
  posterior_overall_effect_distribution <- colMeans(effect_size_posterior, na.rm = T)
  
  # --- 4. Summarize the Overall Average Effect Distribution ---
  
  overall_mean <- mean(posterior_overall_effect_distribution, na.rm = T)
  overall_ci_lower <- quantile(posterior_overall_effect_distribution, probs = 0.025, na.rm = T)
  overall_ci_upper <- quantile(posterior_overall_effect_distribution, probs = 0.975, na.rm = T)
  
  # --- 5. Prepare Data for Plotting ---
  
  plot_data <- tibble(
    locality_id = factor(1:n_localities),
    mean = locality_mean,
    lower_ci = locality_ci_lower,
    upper_ci = locality_ci_upper
  ) |>
    # Order localities by their median effect size for visualization
    arrange(mean) |> 
    mutate(locality_id = fct_inorder(locality_id)) |> 
    slice_sample(n = 10)
  
  # Create a small data frame for the overall effect summary to plot separately
  overall_summary <- tibble(
    mean = overall_mean,
    lower_ci = overall_ci_lower,
    upper_ci = overall_ci_upper
  )
  
  # --- 6. Create Forest Plot ---
  
  forest_plot <- ggplot(plot_data, aes(x = mean, y = locality_id)) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "gray60") +
    geom_point(shape = 15, size = 2, color = "navyblue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_errorbarh(data = overall_summary,
                   aes(y = 0, xmin = lower_ci, xmax = upper_ci), # Place at y=0 or slightly below plot
                   height = 0.4, color = "darkgreen", linewidth = 1) +
    geom_point(data = overall_summary,
               aes(y = 0, x = mean), # Place at y=0 or slightly below plot
               shape = 18, size = 4, color = "darkgreen") +
    # Adjust y-axis limits to make space for the overall summary below the localities
    scale_y_discrete(expand = expansion(add = c(1, 1.5))) +
    labs(
      title = "Forest Plot of Effect Size (Scenario - Baseline )",
      subtitle = paste("Overall effect derived from", n_samples, "posterior samples averaged across", n_localities, "localities"),
      x = paste0("Difference in", metric_nm),
      y = "Sampling Unit ID (Ordered by Mean Effect)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = rel(2), face = "bold"), 
      plot.subtitle = element_text(size = rel(2)),        
      axis.title.x = element_text(size = rel(2)), 
      axis.text.x = element_text(size = rel(2))          
    )
  
  # --- 7. Create gt table ---
  
  table_data_prep <- plot_data |>
    mutate(ci_text = glue("[{sprintf('%.2f', lower_ci)}; {sprintf('%.2f', upper_ci)}]")) |>
    select(locality_id, mean, ci_text)
  
  overall_summary_row <- tibble(
    locality_id = factor("Overall", levels = c(levels(table_data_prep$locality_id), "Overall")),
    mean = overall_mean,
    ci_text = glue("[{sprintf('%.2f', overall_ci_lower)}; {sprintf('%.2f', overall_ci_upper)}]")
  )
  
  # Combine sampling unit and overall summary
  table_data_final <- bind_rows(table_data_prep, overall_summary_row)
  
  summary_table_gt <- table_data_final |>
    gt(rowname_col = "locality_id") |>
    cols_label(
      mean = html("Mean<br>Effect"),
      ci_text = html("95% CrI")
    ) |>
    fmt_number(
      columns = mean,
      decimals = 2
    ) |>
    cols_align(
      align = "center",
      columns = everything()
    ) |>
    tab_style(
      style = list(
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        rows = locality_id == "Overall"
      )
    ) |>
    tab_style(
      style = cell_borders(sides = "top", color = "gray80", weight = px(2)),
      locations = cells_body(rows = locality_id == "Overall")
    ) |>
    # Opciones generales de la tabla
    tab_options(
      table.border.top.width = px(0), # Sin borde superior general
      table.border.bottom.width = px(0), # Sin borde inferior general
      column_labels.border.bottom.width = px(2), # Borde bajo las etiquetas de columna
      column_labels.border.bottom.color = "gray40",
      table_body.border.bottom.color = "#FFFFFF00", # Sin bordes inferiores en el cuerpo
      heading.border.bottom.width = px(0),
      source_notes.border.bottom.width = px(0),
      table.width = pct(100), # Ancho relativo
      data_row.padding = px(12), # Espaciado vertical
      column_labels.padding = px(4)
    ) |>
    opt_table_lines(extent = "none") |> # Empezar sin líneas internas
    opt_horizontal_padding(scale=0.5) # Reducir padding horizontal
  
  
  if(!(metric_nm %in% c("MSA", "SORENSEN","L-PDF*"))){
    #--------------------
    # PPT JYU
    
    localities <- paste0("L", 1:n_samples)
    
    a <-  metric_list[[3]]
    colnames(a) <- localities
    
    a <- a |> 
      as.data.frame() |> 
      pivot_longer(
        cols = everything(),
        names_to = "sample",
        values_to = metric_nm
      ) |> 
      mutate(state = "baseline")
    
    a2 <- metric_list[[4]] 
    colnames(a2) <-  localities
    
    a2 <- a2 |>
      as.data.frame() |> 
      rename_with(.cols = everything(), .fn = ~gsub("V", "L", .x)) |>
      pivot_longer(
        cols = everything(),
        names_to = "sample",
        values_to = metric_nm
      ) |> 
      mutate(state = "scenario")
    
    sample_localities <- sample(localities, 5, replace = F)
    
    p2 <- bind_rows(a, a2) |> 
      filter(sample %in% sample_localities) |> 
      ggplot() + 
      geom_density(mapping = aes(alpha_richness, fill = state), alpha = 0.4) + 
      facet_wrap(~sample, scales = "free", nrow = 5, ncol = 1)  
    
     res[[i]] <- list(metric_nm = metric_nm, overall = overall_summary, forest_plot, table_data_final, summary_table_gt, p2 = p2)
    next()
  }
  
  res[[i]] <- list(metric_nm = metric_nm, overall = overall_summary, forest_plot, table_data_final, summary_table_gt, p2 = NULL)
  
}

names(res) <- metric_nms


