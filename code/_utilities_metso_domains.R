library(tidyverse)
library(modEvA)
library(mop)

# 1. Setup

source(file.path("code", "config_model.r"))

model_nm <- run_config$model_id

fbs_list <- readRDS(
  file.path(
    "data", "covariates",
    paste0("XData_hmsc_coords_", model_nm, ".rds")
  )
)

control_list <- readRDS(
  file.path(
    "data", "covariates",
    paste0("XData_hmsc_control_", model_nm, ".rds")
  )
)

metso_list <- readRDS(
  file.path(
    "data", "covariates",
    paste0("XData_hmsc_metso_", model_nm, ".rds")
  )
)

# 2. Load/extract datasets

fbs <- fbs_list[["XData"]]

control <- control_list[["XData"]]
metso <- metso_list[["XData"]]

treat_control <- bind_rows(control, metso)

projection_group <- c(
  rep("BAU", nrow(control)),
  rep("METSO", nrow(metso))
)

treat_control_poly_ids <- c(
  control_list[["polygon_id"]],
  metso_list[["polygon_id"]]
)

# 3. Define environmental covariates to use for MOP

cov_names <- c(
  "tree_extent",
  "average_stand_diameter",
  "canopy_cover_broadleaves",
  "tree_remove",
  "tree_height_stdev",
  "mean_temperature_yr",
  "stand_age",
  "mean_temperature_yr_last",
  "average_stand_diameter_gt15_mean",
  "average_stand_diameter_gt20_mean",
  "average_stand_length_gt200_mean",
  "canopy_cover_whole_stand_gt60_mean",
  "volume_spruce_gt150"
)

missing_fbs <- setdiff(cov_names, colnames(fbs))
missing_projection <- setdiff(cov_names, colnames(treat_control))

# Subset datasets
ref_data <- fbs |>
  select(all_of(cov_names))

proj_data <- treat_control |>
  select(all_of(cov_names))


# Missing-value check
message(
  "Missing values in reference data: ",
  sum(is.na(ref_data))
)

message(
  "Missing values in projection data: ",
  sum(is.na(proj_data))
)

# 5. MOP calculation

# Scaled Euclidean distance.
# Each projection point is compared with the nearest 10% of
# reference FBS environments.

mop_out <- mop::mop(
  m = as.matrix(ref_data),
  g = as.matrix(proj_data),
  type = "detailed",
  calculate_distance = TRUE,
  where_distance = "all",
  distance = "euclidean",
  scale = TRUE,
  center = TRUE,
  percentage = 10
)

# 6. Overall extrapolation diagnostics

## A. Overall strict extrapolation rate

total_stands <- length(mop_out$mop_basic)

strict_extrap_count <- sum(!is.na(mop_out$mop_basic))

strict_extrap_pct <- (
  strict_extrap_count / total_stands
) * 100

message("Total projected forest stands: ", total_stands)

message(
  "Stands in strict extrapolation: ",
  strict_extrap_count,
  " (",
  round(strict_extrap_pct, 2),
  "%)"
)

# B environmental-distance distribution

dist_summary <- summary(mop_out$mop_distances)

print(dist_summary)


# C. Covariate-specific strict extrapolation

# Extract detailed matrices
low_matrix <- mop_out$mop_detailed$towards_low_end
high_matrix <- mop_out$mop_detailed$towards_high_end

# Convert NAs to zero for rate calculations
low_matrix[is.na(low_matrix)] <- 0
high_matrix[is.na(high_matrix)] <- 0


# Overall variable-specific extrapolation rates
extrap_by_var <- data.frame(
  variable = colnames(low_matrix),
  Low_End_Pct = colMeans(low_matrix > 0) * 100,
  High_End_Pct = colMeans(high_matrix > 0) * 100
) |>
  mutate(
    Total_Univariate_Pct = Low_End_Pct + High_End_Pct
  ) |>
  arrange(desc(Total_Univariate_Pct))

print(extrap_by_var)


# D. Number of simultaneously out-of-range variables


simple_counts <- table(
  mop_out$mop_simple,
  useNA = "ifany"
)

names(simple_counts)[is.na(names(simple_counts))] <- "0 (In Range)"

print(simple_counts)


# E. Construct stand-level MOP diagnostic dataset

mop_diagnostics <- tibble(
  stand_id = treat_control_poly_ids,
  group = factor(
    projection_group,
    levels = c("BAU", "METSO")
  ),
  
  # 1 = at least one predictor outside FBS reference range
  mop_strict_extrap = ifelse(
    is.na(mop_out$mop_basic),
    0L,
    1L
  ),
  
  # Continuous multivariate environmental distance
  mop_distance = mop_out$mop_distances,
  
  # Number of predictors outside FBS range
  number_out_vars = ifelse(
    is.na(mop_out$mop_simple),
    0L,
    mop_out$mop_simple
  )
)

# 7. Compare MOP diagnostics between BAU and METSO

mop_summary_by_group <- mop_diagnostics |>
  group_by(group) |>
  summarise(
    n = n(),
    n_strict_extrap = sum(mop_strict_extrap),
    strict_extrap_pct = mean(mop_strict_extrap) * 100,
    mean_distance = mean(
      mop_distance,
      na.rm = TRUE
    ),
    sd_distance = sd(
      mop_distance,
      na.rm = TRUE
    ),
    q25_distance = quantile(
      mop_distance,
      probs = 0.25,
      na.rm = TRUE
    ),
    median_distance = median(
      mop_distance,
      na.rm = TRUE
    ),
    q75_distance = quantile(
      mop_distance,
      probs = 0.75,
      na.rm = TRUE
    ),
    q95_distance = quantile(
      mop_distance,
      probs = 0.95,
      na.rm = TRUE
    ),
    max_distance = max(
      mop_distance,
      na.rm = TRUE
    ),
    mean_number_out_vars = mean(
      number_out_vars,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(mop_summary_by_group)

# B. Distribution of number of out-of-range variables
#     separately for BAU and METSO

out_vars_by_group <- mop_diagnostics |>
  count(
    group,
    number_out_vars,
    name = "n"
  ) |>
  group_by(group) |>
  mutate(
    pct = n / sum(n) * 100
  ) |>
  ungroup()

print(out_vars_by_group)

# C. Covariate-specific extrapolation by BAU/METSO group

# Function to calculate variable-specific extrapolation
# within one projection group
calculate_extrap_by_group <- function(group_name) {
  
  idx <- projection_group == group_name
  
  tibble(
    group = group_name,
    variable = colnames(low_matrix),
    
    Low_End_Pct =
      colMeans(
        low_matrix[idx, , drop = FALSE] > 0
      ) * 100,
    
    High_End_Pct =
      colMeans(
        high_matrix[idx, , drop = FALSE] > 0
      ) * 100
  ) |>
    mutate(
      Total_Univariate_Pct =
        Low_End_Pct + High_End_Pct
    )
}

extrap_by_group_var <- bind_rows(
  calculate_extrap_by_group("BAU"),
  calculate_extrap_by_group("METSO")
) |>
  arrange(
    variable,
    group
  )

print(extrap_by_group_var)


# D. Direct BAU vs METSO comparison for each covariate

extrap_group_comparison <- extrap_by_group_var |>
  select(
    group,
    variable,
    Total_Univariate_Pct
  ) |>
  pivot_wider(
    names_from = group,
    values_from = Total_Univariate_Pct
  ) |>
  mutate(
    METSO_minus_BAU_pp = METSO - BAU
  ) |>
  arrange(
    desc(abs(METSO_minus_BAU_pp))
  )

# E. Compare low-end vs high-end extrapolation directly

extrap_group_detailed <- extrap_by_group_var |>
  select(group, variable, Low_End_Pct, High_End_Pct, Total_Univariate_Pct
  ) |>
  pivot_wider(
    names_from = group,
    values_from = c(Low_End_Pct, High_End_Pct, Total_Univariate_Pct
    )
  ) |>
  mutate(
    Total_METSO_minus_BAU_pp = Total_Univariate_Pct_METSO - Total_Univariate_Pct_BAU,
    Low_METSO_minus_BAU_pp = Low_End_Pct_METSO - Low_End_Pct_BAU,
    High_METSO_minus_BAU_pp = High_End_Pct_METSO -High_End_Pct_BAU
  ) |>
  arrange(desc(abs(Total_METSO_minus_BAU_pp)))

# 8. Descriptive comparison of MOP distance

distance_comparison <- mop_diagnostics |>
  group_by(group) |>
  summarise(
    n = n(),
    mean = mean(
      mop_distance,
      na.rm = TRUE
    ),
    sd = sd(
      mop_distance,
      na.rm = TRUE
    ),
    median = median(
      mop_distance,
      na.rm = TRUE
    ),
    IQR = IQR(
      mop_distance,
      na.rm = TRUE
    ),
    q05 = quantile(
      mop_distance,
      0.05,
      na.rm = TRUE
    ),
    q25 = quantile(
      mop_distance,
      0.25,
      na.rm = TRUE
    ),
    q75 = quantile(
      mop_distance,
      0.75,
      na.rm = TRUE
    ),
    q95 = quantile(
      mop_distance,
      0.95,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Save results

dir.create("results/mop", showWarnings = FALSE, recursive = TRUE)

# Stand-level diagnostic file
write.csv(mop_diagnostics, "results/mop/mop_diagnostics_mask.csv", row.names = FALSE)

# Overall covariate-specific extrapolation
write.csv(extrap_by_var, "results/mop/mop_extrapolation_by_covariate.csv", row.names = FALSE)

# Group-level MOP summary
write.csv(mop_summary_by_group, "results/mop/mop_summary_by_group.csv", row.names = FALSE)

# Number of out-of-range variables by group
write.csv(out_vars_by_group, "results/mop/mop_number_out_vars_by_group.csv", row.names = FALSE)

# Covariate-specific extrapolation separately for BAU/METSO
write.csv(extrap_by_group_var, "results/mop/mop_extrapolation_by_covariate_group.csv", row.names = FALSE)

# Direct METSO-vs-BAU differences
write.csv(extrap_group_detailed, "results/mop/mop_extrapolation_METSO_vs_BAU.csv",row.names = FALSE)

# Continuous MOP-distance comparison
write.csv(distance_comparison, "results/mop/mop_distance_by_group.csv", row.names = FALSE)