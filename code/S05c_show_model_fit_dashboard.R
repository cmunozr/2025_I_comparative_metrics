# J-SDM Model Diagnostic Builder & Interactive Dashboard

library(Hmsc)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(shiny)
library(ggplot2)
library(plotly)

# --- 1. Paths & Configurations ---
base_dir          <- "models"
species_info_path <- "data/species_info_table.csv"
master_fit_path   <- file.path(base_dir, "model_fit_master_all_metrics.csv")

pa_model_dir  <- file.path(base_dir, "fbs_M016PA_thin_150_samples_1000_chains_4")
acp_model_dir <- file.path(base_dir, "fbs_M016_thin_250_samples_1000_chains_4")


# Target focal model identifier
focal_model_id <- "fbs_M016"

# Mapping strategy identifiers to clean column suffixes
strategy_code_map <- c(
  "Full_Model"       = "fit",
  "route_blocked_cv" = "cv_route",
  "metso_holdout"    = "ho_metso",
  "north_south"      = "ho_north"
)

# --- 2. Load Base Species Data ---
table_info_species <- read.csv(species_info_path)
spp <- table_info_species$species

# --- 3. Extract All Strategies from Master CSV ---
if (!file.exists(master_fit_path)) {
  stop("Master fit CSV not found. Run S05b_show_model_fit_grouped.R first.")
}

master_fit <- read_csv(master_fit_path, show_col_types = FALSE)

# Standardize strategy labels
master_fit_clean <- master_fit %>%
  filter(Model == focal_model_id) %>%
  mutate(Strat_Code = recode(as.character(Strategy), !!!strategy_code_map, .default = as.character(Strategy)))

# A. Presence-Absence (PA) across all strategies
df_PA_all_strat <- master_fit_clean %>%
  filter(Model_Type == "Presence-Absence") %>%
  select(Species, Strat_Code, any_of(c("AUC", "TjurR2", "RMSE"))) %>%
  pivot_wider(
    names_from = Strat_Code,
    values_from = any_of(c("AUC", "TjurR2", "RMSE")),
    names_glue = "PA_{.value}_{Strat_Code}"
  ) %>%
  rename(species = Species)

# B. Conditional Abundance (aCp) across all strategies
df_aCp_all_strat <- master_fit_clean %>%
  filter(Model_Type == "Continuous Abundance") %>%
  select(Species, Strat_Code, any_of(c("SR2", "RMSE"))) %>%
  pivot_wider(
    names_from = Strat_Code,
    values_from = any_of(c("SR2", "RMSE")),
    names_glue = "aCp_{.value}_{Strat_Code}"
  ) %>%
  rename(species = Species)

# --- 5. Helper Function: Extract Variance Partitioning & Betas ---
extract_parameter_estimates <- function(model_dir, suffix) {
  vp_file   <- file.path(model_dir, "parameter_estimates", "parameter_estimates_VP.csv")
  beta_file <- file.path(model_dir, "parameter_estimates", "parameter_estimates_Beta.xlsx")
  
  df_vp <- NULL
  if (file.exists(vp_file)) {
    vp_raw <- read.csv(vp_file, check.names = FALSE, row.names = 1)
    df_vp <- as.data.frame(t(vp_raw))
    
    fixed_cols <- intersect(c("climate", "forest_structure", "survey"), colnames(df_vp))
    df_vp[[paste0("Fixed_", suffix)]] <- rowSums(df_vp[, fixed_cols, drop = FALSE])
    
    df_vp <- df_vp %>%
      select(-any_of(c("Intercept", "climate", "forest_structure", "survey"))) %>%
      rename_with(~paste0(.x, "_", suffix), .cols = everything())
    
    df_vp$species <- rownames(df_vp)
  }
  
  df_beta <- NULL
  if (file.exists(beta_file)) {
    df_beta <- read.xlsx(beta_file) %>%
      mutate(
        beta_95 = if_else(support >= 0.95, 1, NA_real_),
        beta_85 = if_else(support >= 0.85, 1, NA_real_)
      ) %>%
      rename(species = Species) %>%
      filter(!is.na(beta_85)) %>%
      group_by(species) %>%
      summarize(
        !!paste0("sum_beta_95_", suffix) := sum(beta_95, na.rm = TRUE),
        !!paste0("sum_beta_85_", suffix) := sum(beta_85, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  if (!is.null(df_vp) && !is.null(df_beta)) {
    return(full_join(df_vp, df_beta, by = "species"))
  } else if (!is.null(df_vp)) {
    return(df_vp)
  } else {
    return(df_beta)
  }
}

# --- 6. Extract Parameters for Both Sub-Models ---
params_PA  <- extract_parameter_estimates(pa_model_dir, suffix = "PA")
params_aCp <- extract_parameter_estimates(acp_model_dir, suffix = "aCp")

# --- 7. Assemble Master Diagnostic Table ---
diagnostic_df <- table_info_species %>%
  left_join(df_PA_all_strat, by = "species") %>%
  left_join(df_aCp_all_strat, by = "species")

if (!is.null(params_PA)) {
  diagnostic_df <- diagnostic_df %>% left_join(params_PA, by = "species")
}

if (!is.null(params_aCp)) {
  diagnostic_df <- diagnostic_df %>% left_join(params_aCp, by = "species")
}

diagnostic_df <- diagnostic_df %>%
  mutate(across(starts_with("sum_beta"), ~tidyr::replace_na(.x, 0)))

# Save comprehensive table
output_csv <- file.path(base_dir, "diagnostic_models_df.csv")
write.csv(diagnostic_df, output_csv, row.names = FALSE)
message("Saved multi-strategy diagnostic dataset to: ", output_csv)

# =========================================================================
# 8. Interactive Shiny Diagnostic Application
# =========================================================================

df <- read.csv(output_csv)

numeric_cols     <- names(df)[sapply(df, is.numeric)]
categorical_cols <- c("None", names(df)[sapply(df, is.character)])

ui <- fluidPage(
  titlePanel("J-SDM Multi-Strategy Diagnostic Dashboard"),
  sidebarLayout(
    sidebarPanel(
      h4("Axis Selection"),
      selectInput("x_var", "X-Axis Variable:", 
                  choices = numeric_cols, 
                  selected = if ("PA_AUC_cv_route" %in% numeric_cols) "PA_AUC_cv_route" else numeric_cols[1]),
      
      selectInput("y_var", "Y-Axis Variable:", 
                  choices = numeric_cols, 
                  selected = if ("PA_AUC_ho_metso" %in% numeric_cols) "PA_AUC_ho_metso" else numeric_cols[2]),
      
      h4("Visual Grouping"),
      selectInput("color_var", "Color by Trait / Classification:", 
                  choices = categorical_cols, 
                  selected = if ("Specialist" %in% categorical_cols) "Specialist" else "None"),
      
      h4("Data Filtering"),
      sliderInput("prev_filter", "Minimum Prevalence (%):", min = 0, max = 100, value = 0, step = 1),
      hr(),
      helpText("Evaluate transferability: Compare explanatory fits, spatial cross-validation (cv_route), and METSO holdouts (ho_metso).")
    ),
    mainPanel(
      plotlyOutput("scatterPlot", height = "600px"),
      br(),
      h4("Summary Statistics by Group"),
      tableOutput("summaryStats"),
      br(),
      dataTableOutput("dataTable")
    )
  )
)

server <- function(input, output) {
  filtered_data <- reactive({
    if ("prev_clean" %in% colnames(df)) {
      df %>% filter(prev_clean >= input$prev_filter)
    } else {
      df
    }
  })
  
  output$scatterPlot <- renderPlotly({
    req(input$x_var, input$y_var)
    p <- ggplot(filtered_data(), aes_string(x = input$x_var, y = input$y_var, text = "species")) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("Validation Comparison:", input$x_var, "vs", input$y_var),
        x = input$x_var,
        y = input$y_var
      )
    
    if (input$color_var != "None") {
      p <- p + aes_string(color = input$color_var) +
        scale_color_viridis_d(option = "plasma", na.value = "grey50")
    }
    
    # 1:1 Reference line for comparing validation holdouts
    if ((grepl("AUC", input$x_var) && grepl("AUC", input$y_var)) || 
        (grepl("SR2", input$x_var) && grepl("SR2", input$y_var)) ||
        (grepl("RMSE", input$x_var) && grepl("RMSE", input$y_var)) ||
        (grepl("TjurR2", input$x_var) && grepl("TjurR2", input$y_var))) {
      p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
    }
    
    ggplotly(p, tooltip = "text")
  })
  
  output$summaryStats <- renderTable({
    req(input$x_var, input$y_var, input$color_var)
    data <- filtered_data()
    
    calc_stats <- function(df_in, var_name, group_col = NULL) {
      if (is.null(group_col)) {
        df_in %>%
          summarise(
            Group    = "All Species (Total)",
            Variable = var_name,
            Mean     = mean(.data[[var_name]], na.rm = TRUE),
            SD       = sd(.data[[var_name]], na.rm = TRUE),
            Min      = min(.data[[var_name]], na.rm = TRUE),
            Max      = max(.data[[var_name]], na.rm = TRUE)
          )
      } else {
        df_in %>%
          group_by(.data[[group_col]]) %>%
          summarise(
            Variable = var_name,
            Mean     = mean(.data[[var_name]], na.rm = TRUE),
            SD       = sd(.data[[var_name]], na.rm = TRUE),
            Min      = min(.data[[var_name]], na.rm = TRUE),
            Max      = max(.data[[var_name]], na.rm = TRUE),
            .groups  = "drop"
          ) %>%
          rename(Group = !!sym(group_col)) %>%
          mutate(Group = as.character(Group))
      }
    }
    
    total_x <- calc_stats(data, input$x_var)
    total_y <- calc_stats(data, input$y_var)
    summary_table <- bind_rows(total_x, total_y)
    
    if (input$color_var != "None") {
      grouped_x <- calc_stats(data, input$x_var, input$color_var)
      grouped_y <- calc_stats(data, input$y_var, input$color_var)
      summary_table <- bind_rows(summary_table, grouped_x, grouped_y)
    }
    
    summary_table %>% arrange(Variable, Group)
  }, digits = 4)
  
  output$dataTable <- renderDataTable({
    filtered_data() %>%
      select(species, input$x_var, input$y_var, any_of(input$color_var))
  }, options = list(pageLength = 10))
}

shinyApp(ui = ui, server = server)