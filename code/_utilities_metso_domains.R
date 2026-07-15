library(tidyverse)
library(modEvA)
library(mop)

# 1. Source configuration and load datasets
source(file.path("code", "config_model.r"))

model_nm <- run_config$model_id 

# read complete object
fbs_list <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_coords_", model_nm, ".rds")))
control_list <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_control_", model_nm, ".rds")))
metso_list <- readRDS(file.path("data", "covariates", paste0("XData_hmsc_metso_", model_nm, ".rds")))

# extract data
fbs <- fbs_list[["XData"]]
treat_control <- bind_rows(control_list[["XData"]], metso_list[["XData"]])

# extract id
treat_control_poly_ids <- c(control_list[["polygon_id"]], metso_list[["polygon_id"]]) |> 
  droplevels()

# 2. Dynamically identify model covariates

cov_names <- intersect(colnames(fbs), colnames(treat_control))

message("Running diagnostics on ", length(cov_names), " variables")

# Subset datasets
ref_data <- fbs |>  
  select(all_of(cov_names))
proj_data <- treat_control |> 
  select(all_of(cov_names))

# 3. Vectorized MESS Calculation Function

mess_out <- modEvA::MESS(
  V = as.data.frame(ref_data), 
  P = as.data.frame(proj_data))

# Step 4. MOP Calculation 

mop_out <- mop::mop(
  m = as.matrix(ref_data),
  g = as.matrix(proj_data),
  type = "detailed",           # "basic" calculates strict extrapolation and distance
  calculate_distance = TRUE,   # Quantifies the environmental distance to the envelope
  where_distance = "all",      # Computes distance values for all projection points
  distance = "mahalanobis",    # Distance metric
  scale = TRUE,                # Standardizes variables (essential for different units)
  center = TRUE,               # Centers variables
  percentage = 10,             # compare the projection points to the closest 10% of reference environments
)

# Step 5: Analyze Diagnostics

## A. Overall Extrapolation Rates
total_stands <- length(mop_out$mop_basic)
strict_extrap_count <- sum(!is.na(mop_out$mop_basic))
strict_extrap_pct <- (strict_extrap_count / total_stands) * 100

message("Total projected forest stands: ", total_stands)
message("Stands in strict extrapolation: ", strict_extrap_count, " (", round(strict_extrap_pct, 2), "%)")

## B. Environmental Distance Distribution, multivariate distance (perypherality) statistics
dist_summary <- summary(mop_out$mop_distances)
print(dist_summary)

# Identify a threshold for high peripherality (e.g., 95th percentile)
distance_95pct <- quantile(mop_out$mop_distances, 0.95, na.rm = TRUE)
message("95th percentile distance threshold: ", round(distance_95pct, 4))

## C. Deconstruct Covariate-Specific Drivers (Low vs. High End)
# Extract detailed matrices
low_matrix <- mop_out$mop_detailed$towards_low_end
high_matrix <- mop_out$mop_detailed$towards_high_end

# Convert NAs to 0 for column-wise rate calculations
low_matrix[is.na(low_matrix)] <- 0
high_matrix[is.na(high_matrix)] <- 0

# Compile variable-specific extrapolation rates
extrap_by_var <- data.frame(
  Low_End_Pct = colMeans(low_matrix > 0) * 100,
  High_End_Pct = colMeans(high_matrix > 0) * 100
)

# Calculate total univariate extrapolation rate per variable
extrap_by_var$Total_Univariate_Pct <- extrap_by_var$Low_End_Pct + extrap_by_var$High_End_Pct
extrap_by_var <- extrap_by_var[order(-extrap_by_var$Total_Univariate_Pct), ]

# Extrapolation rate by individual covariate
print(extrap_by_var)

## D. Variable Intersections
# Extract counts of how many variables fail simultaneously per stand
simple_counts <- table(mop_out$mop_simple, useNA = "ifany")
names(simple_counts)[is.na(names(simple_counts))] <- "0 (In Range)"

#Number of out-of-range variables per stand
print(simple_counts)

# E. Export Diagnostic Mask for HMSC Prediction Scripts
# Link to your projection identifiers (e.g., stand_id)
mop_diagnostics <- data.frame(
  stand_id = treat_control_poly_ids,
  mop_strict_extrap = ifelse(is.na(mop_out$mop_basic), 0, 1),
  mop_distance = mop_out$mop_distances,
  number_out_vars = ifelse(is.na(mop_out$mop_simple), 0, mop_out$mop_simple)
)

# Flag highly peripheral stands (inside range but in the top 5% of environmental distance)
mop_diagnostics$is_highly_peripheral <- ifelse(
  mop_diagnostics$mop_distance >= distance_95pct & mop_diagnostics$mop_strict_extrap == 0, 
  1, 
  0
)

write.csv(mop_diagnostics, "data/metso/mop_diagnostics_mask.csv", row.names = FALSE)

