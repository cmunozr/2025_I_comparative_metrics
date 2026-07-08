library(data.table)
library(dplyr)
library(arrow)

# Load data assets
load(file.path("data", "metso", "raw2", "treatment_long.RData"))
matched_units <- arrow::read_parquet("data/metso/raw2/matched_units.parquet")
sp <- arrow::read_parquet("data/metso/raw2/stacked_panel.parquet")

# 1. Data description
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

# 2: Process the stacked panel data
sp <- data.table::as.data.table(sp)

P <- 5L # Pre-event window years
L <- 5L # Post-event window snapshot (5-year target horizon)
first_data_year <- 2001
last_data_year <- 2021 

# Extract cohort enrollment years within the data.table class
cohort_years <- sp[
  treated_stack == 1L,
  .(cohort_year = as.integer(min(year, na.rm = TRUE))),
  by = stack_id
]

# Identify cohorts with complete historical and post-treatment records
eligible_cohorts <- cohort_years[
  cohort_year - P >= first_data_year & 
  cohort_year + L <= last_data_year
]

# Append the explicit calendar year corresponding to event_time == 5
eligible_cohorts <- eligible_cohorts[, year_event_5 := cohort_year + L] |> 
  dplyr::filter(year_event_5 %in% c(2009, 2011, 2013, 2015, 2017, 2019, 2021))

# 3: Filter matched units and join temporal attributes
matched_units <- data.table::as.data.table(matched_units)

# Filter for eligible stacks
matched_units_filtered <- matched_units[
  stack_id %in% eligible_cohorts$stack_id
]

# Join the target event year back to both treated and control units
matched_units_filtered <- merge(
    matched_units_filtered,
    eligible_cohorts[, .(stack_id, year_event_5)],
    by = "stack_id",
    all.x = TRUE
  )|> 
    dplyr::mutate(
      metso = if_else(match_role == "treated", 1, 0)  
    )

data.table::fwrite(
  matched_units_filtered, 
  file.path("data", "metso", "treatment_control_standid.csv")
)

