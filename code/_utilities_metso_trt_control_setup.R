library(data.table)
library(dplyr)
library(arrow)

# Load data assets
load(file.path("data", "metso", "raw2", "treatment_long.RData"))
matched_units <- arrow::read_parquet("data/metso/raw2/matched_units.parquet")
sp <- arrow::read_parquet("data/metso/raw2/stacked_panel.parquet")

# 0. Data description
trt_stands_summary <- trt.long |> 
  dplyr::filter(metso_any_bin != 0, year <= 2021) |> 
  dplyr::mutate(across(
    starts_with("metso_") & !contains("any"), 
    ~ if_else(. > 0, 1, 0)
  )) |> 
  dplyr::mutate(contract_type = case_when(
    metso_perm == 1       ~ "metso_perm",
    metso_statepurch == 1 ~ "metso_statepurch",
    metso_temp == 1       ~ "metso_temp",
    metso_10k_13k == 1    ~ "metso_10k_13k",
    metso_kemera == 1     ~  "metso_kemera",
    TRUE                  ~ NA_character_
  )) |> 
  dplyr::group_by(standid) |> 
  dplyr::summarise(
    year_in = min(year), 
    contract_type = unique(na.omit(contract_type))[1],
    .groups = "drop"
  ) |> 
  dplyr::mutate(year_length = 2021 - year_in, metso = 1)

table(trt_stands_summary$year_in) |>barplot()

table(trt_stands_summary$contract_type) |> barplot()

#--------------

# 1: Process the stacked panel data
sp <- data.table::as.data.table(sp)

P <- 1L # Pre-event window years
first_data_year <- 2009
last_data_year  <- 2021

# Extract cohort enrollment years within the data.table class
cohort_years <- sp[
  treated_stack == 1L,
  .(stack_id, cohort_year = unique(g)),
  by = stack_id
]

# Process 5-year post-event horizons
eligible_5 <- cohort_years[
  cohort_year - P >= first_data_year &
    cohort_year + 5 <= last_data_year
]

eligible_5[, `:=`(event_horizon = 5L, year_event = cohort_year + 5L)]

# Process 6-year post-event horizons
eligible_6 <- cohort_years[
  cohort_year - P >= first_data_year &
    cohort_year + 6 <= last_data_year
]

eligible_6[, `:=`(event_horizon = 6L, year_event = cohort_year + 6L)]

# Combine both horizons into a unified tracking data table
eligible_cohorts <- rbind(eligible_5, eligible_6) 

# 2: Filter matched units and join temporal attributes
matched_units <- data.table::as.data.table(matched_units)

# Merge mapping metrics to duplicate stacks qualifying for both targets
matched_units_filtered <- matched_units |> 
  merge(
    eligible_cohorts[, .(stack_id, event_horizon, year_event)],
    by = "stack_id",
    all.y = TRUE
  )|> 
  dplyr::filter(
    stack_id %in% eligible_cohorts$stack_id
  ) 

# Allocate explicit treatment indicators using data.table infrastructure
matched_units_filtered[, metso := data.table::fifelse(match_role == "treated", 1L, 0L)]

matched_units_filtered <- matched_units_filtered |> 
  as.data.frame() |> 
  dplyr::rename(year = year_event) |>
  dplyr::mutate(
    standid = paste0(standid, "_", year, "_", event_horizon)
  ) |> 
  dplyr::filter(
    year %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021)
  )

data.table::fwrite(
  matched_units_filtered, 
  file.path("data", "metso", "treatment_control_standid.csv")
)
