#' DESCRIPTION:
#' Creates a descriptive table comparing covariates in METSO and BAU stands
#' for the selected hurdle model.
#'
#' METHODOLOGY:
#' - Selects model covariates marked for the summary table in the dictionary.
#' - Summarises continuous variables as means and standard deviations and
#'   binary variables as counts and percentages.
#' - Reports METSO minus BAU differences and standardised mean differences.
#'
#' INPUTS:
#' - Model selection: code/config_model.R
#' - Covariate keys: models/[selected_model]/parameter_estimates/covariate_key_table.csv
#' - Covariate dictionary: data/covariates/mapping_covariate_functions.xlsx
#' - Stand covariates: data/covariates/XData_hmsc_[metso/control]_[model_id].rds
#' - Survey covariates: data/covariates/XData_hmsc_coords_[model_id].rds
#'
#' OUTPUTS:
#' - data/covariates/summary_table.html
#' - data/covariates/summary_table.tex

library(dplyr)
library(gtsummary)
library(gt)
library(openxlsx)

source(file.path("code", "config_model.R"))

#-----------
# Raw METSO - BAU difference for continuous variables
raw_diff_continuous <- function(data, variable, by, ...) {
  diff <- mean(
    data[[variable]][data[[by]] == "METSO"],
    na.rm = TRUE
  ) -
    mean(
      data[[variable]][data[[by]] == "BAU"],
      na.rm = TRUE
    )
  tibble::tibble(raw_diff = diff)
}

# Raw METSO - BAU difference for binary variables
# Returned in percentage points
raw_diff_binary <- function(data, variable, by, ...) {
  p_metso <- mean(
    data[[variable]][data[[by]] == "METSO"] == "Yes",
    na.rm = TRUE
  )
  p_bau <- mean(
    data[[variable]][data[[by]] == "BAU"] == "Yes",
    na.rm = TRUE
  )
  tibble::tibble(
    raw_diff = 100 * (p_metso - p_bau)
  )
}
#---------------------

model_nm <- run_config$model_id 

# 1. Load Dictionary and Identify Variables
dict_covar <- read.xlsx("data/covariates/mapping_covariate_functions.xlsx")

# Load variables used in models
pa_cov_keys <- file.path("models", selected_models$PA, "parameter_estimates", "covariate_key_table.csv") |> 
  read.csv()
ab_cov_keys <- file.path("models", selected_models$aCp, "parameter_estimates", "covariate_key_table.csv") |> 
  read.csv()

cov_keys <- bind_rows(pa_cov_keys, ab_cov_keys) |> 
  filter(Domain != "Intercept")
cov_keys <- cov_keys[!duplicated(cov_keys),]

# Filter dictionary to variables included in model
dict_covar <- dict_covar |> 
  filter(var_target %in% cov_keys$Model_Covariate)

# Include both continuous and binary target variables
is_for_summ_table <- dict_covar$var_target[dict_covar$to_table_summary == 1]

# Identify binary variables explicitly from dictionary
is_binary <- dict_covar$var_target[dict_covar$to_table_summary == 1 & dict_covar$continuos == 0]

# Build the named label vector for tbl_summary
# Clean potential unit formatting artifacts from Excel
dict_labels <- dict_covar |>
  filter(to_table_summary == 1) |>
  mutate(
    clean_label = Inline_Abbreviated |>
      gsub(pattern = "m_3/ha", replacement = "m³/ha", fixed = TRUE) |>
      gsub(pattern = "\\%", replacement = "%", fixed = TRUE)
  )

var_label_list <- setNames(as.list(dict_labels$clean_label), dict_labels$var_target)

# 2. Read Model Covariates
control_list <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_control_", model_nm, ".rds")))
metso_list   <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_metso_", model_nm, ".rds")))

for (group in c("control", "metso")) {
  group_data <- if (group == "control") control_list$XData else metso_list$XData
  missing_vars <- setdiff(is_for_summ_table, names(group_data))
  if (length(missing_vars) > 0) {
    stop("Missing table covariates in ", group, ": ", paste(missing_vars, collapse = ", "))
  }
}

treated <- metso_list[["XData"]] |> 
  select(all_of(is_for_summ_table))

control <- control_list[["XData"]] |> 
  select(all_of(is_for_summ_table))

# 3. Combine Datasets and Format Factors
combined_data <- bind_rows(
  control |> mutate(Group = "BAU"),
  treated |> mutate(Group = "METSO")
) |>
  mutate(
    Group = factor(Group, levels = c("METSO", "BAU")),
    across(
      any_of(is_binary),
      ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))
    )
  )

# 4. Build gtsummary Table
formatted_table <- combined_data |>
  tbl_summary(
    by = Group,
    label = var_label_list,
    type = list(
      all_continuous() ~ "continuous",
      any_of(is_binary) ~ "dichotomous"
    ),
    value = list(
      any_of(is_binary) ~ "Yes"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ c(3, 3),
      all_categorical() ~ c(0, 1)
    )
  ) |>
  add_stat(
    fns = list(
      all_continuous() ~ raw_diff_continuous,
      any_of(is_binary) ~ raw_diff_binary
    )
  ) |>
  add_difference(
    test = everything() ~ "smd"
  ) |>
  modify_table_body(
    ~ .x |> 
      dplyr::select(-conf.low, -conf.high)
  ) |>
  modify_spanning_header(
    all_stat_cols() ~ "**5–6 years post-METSO enrollment**"
  ) |>
  modify_header(
    label ~ "**Covariate**",
    raw_diff ~ "**METSO − BAU**",
    estimate ~ "**SMD**"
  ) |>
  modify_fmt_fun(
    raw_diff ~ label_style_number(digits = 2),
    estimate ~ label_style_number(digits = 2)
  ) |> 
  modify_footnote(
    all_stat_cols() ~ 
      "METSO − BAU is the difference in means for continuous variables and the difference in prevalence (percentage points) for binary variables."
  )

# 4b. Add FBS route-year summaries
fbs_list <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_coords_", model_nm, ".rds")))
missing_fbs_vars <- setdiff(is_for_summ_table, names(fbs_list$XData))
if (length(missing_fbs_vars) > 0) {
  stop("Missing table covariates in FBS: ", paste(missing_fbs_vars, collapse = ", "))
}

fbs_data <- fbs_list$XData |>
  select(all_of(is_for_summ_table)) |>
  mutate(across(any_of(is_binary), ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))))

fbs_summary <- fbs_data |>
  tbl_summary(
    type = list(all_continuous() ~ "continuous", any_of(is_binary) ~ "dichotomous"),
    value = list(any_of(is_binary) ~ "Yes"),
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ c(3, 3), all_categorical() ~ c(0, 1))
  )

formatted_table <- formatted_table |>
  modify_table_body(
    ~ .x |>
      left_join(select(fbs_summary$table_body, variable, row_type, stat_0), by = c("variable", "row_type")) |>
      relocate(stat_0, .before = stat_1)
  ) |>
  modify_header(stat_0 ~ paste0("**FBS**  \nN = ", format(nrow(fbs_data), big.mark = ","))) |>
  modify_column_unhide(columns = stat_0)

# 5. Export to HTML and LaTeX (for Overleaf / elsarticle)
gt_table <- formatted_table |> 
  as_gt() |>
  tab_options(
    table.font.size = "14px",
    heading.align = "left"
  )

# HTML preview
gtsave(gt_table, "data/covariates/summary_table.html")

# LaTeX file for direct \input{summary_table.tex} in Overleaf
gtsave(gt_table, "data/covariates/summary_table.tex")
