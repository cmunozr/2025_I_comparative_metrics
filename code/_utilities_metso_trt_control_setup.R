library(data.table)
library(dplyr)

load(file.path("data", "metso", "raw", "treatment_long.RData"))
matched_pairs <- readRDS("data/metso/raw/matched_pairs.rds") # got 2026-01-12

#------------

# data from 2000 to 2024 (25 years)
num_stands <- unique(trt.long$standid) |> 
  length()
# 57285 unique number stands
num_stands * length(unique(trt.long$year))
# number stands * number of years = 1432125

trt.long.filter <- trt.long |> 
  filter(metso_any_bin != 0, year <= 2021)
table(trt.long.filter$year) |> 
  barplot(); abline(h = length(unique(trt.long.filter$standid)), col = "red", lty = 3)

trt.stands <- trt.long.filter |> 
  mutate(across(
    starts_with("metso_") & !contains("any"), 
    ~ if_else(. > 0, 1, 0)
  )) |> 
  mutate(contract_type = case_when(
    metso_perm == 1       ~ "metso_perm",
    metso_statepurch == 1 ~ "metso_statepurch",
    metso_temp == 1       ~ "metso_temp",
    metso_10k_13k == 1    ~ "metso_10k_13k",
    TRUE                  ~ NA_character_ # For rows where all are 0
  )) |> 
  group_by(standid) |> 
  summarise(year_in = min(year), contract_type = unique(na.omit(contract_type))) |> 
  mutate(year_length = 2021-year_in, metso = 1)
table(trt.stands$year_in) |> 
  barplot()
table(trt.stands$contract_type) |> 
  barplot()

matched_pairs <- matched_pairs |> 
  filter(standid_treated %in% trt.stands$standid)

#------------

control.stands <- matched_pairs |>  
  select(standid_matched_control) |> 
  rename(standid = standid_matched_control) |> 
  mutate(metso = 0)

trt.control <- bind_rows(trt.stands, control.stands)

data.table::fwrite(trt.control, file.path("data","metso", "treatment_control_standid.csv"))


