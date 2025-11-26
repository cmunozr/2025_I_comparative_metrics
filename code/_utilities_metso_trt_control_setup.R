library(data.table)
library(dplyr)

load(file.path("data", "metso", "raw", "treatment_long.RData"))

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
#------------

ard <- readRDS(file.path("data", "metso", "raw", "ard_long.RDS"))
# 273427 unique number stands
# number stands * number of years = 5741967

index <- !(ard$standid %in% trt.stands$standid)
control.stands <- ard[index, ] |>  
  filter(year == 2021) |> 
  select(standid) |> 
  mutate(metso = 0)
# it doesnt matter which year, they are the same in each year
# 228210 unique number stands

#------------
# stand_dt_full <- fread("data/metso/stand_dt_full.csv")
# 
# head(stand_dt_full)
# names(stand_dt_full)
# 
# stand_dt_full <- stand_dt_full |> select(dplyr::contains("metso"))
# 
# unique(stand_dt_full$standid) |>  length()
# rm(stand_dt_full);gc()
# 
# # all the stands? 12385834
# # dataset information is covariate

trt.control <- bind_rows(trt.stands, control.stands)

data.table::fwrite(trt.control, file.path("data","metso", "treatment_control_standid.csv"))
