library(tidyverse)
library(Hmsc)
library(parallel)
library(here)
library(arrow)

# setup

os <- Sys.info()['sysname']
cores <- if (os == "Windows") 2 else 10

if(!(Sys.getenv("RSTUDIO") == "1")){
  setwd(here::here())
}

source(file.path("code", "diversity_metrics_functions.R"))
source(file.path("code", "config_model.R"))
name_model <- generate_run_name(run_config)


# 1. Data loading

start_block <- Sys.time()

load(file.path("models", paste0("unfitted_", run_config$model_id, ".RData")))
name_spp <- colnames(models$fbs_M008$Y)

output_dir <- file.path("results", "metrics")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

set.seed(11072024)

mfeval <- readRDS(file.path("models", name_model, "model_fit", paste0("mfeval_", name_model, "_ho1.rds"))) |> 
  pluck("RMSE")
names(mfeval) <- name_spp

TrData <- readRDS(file.path("data", "traits", "TrData_complete.rds"))
rownames(TrData) <- name_spp


# 2. Batch pre-processing


expected_val <- TRUE 
expected_string <- if(expected_val) "_expected_true" else "_expected_false"

## Open arrow dataset and filter

predY_ds <- open_dataset(
  file.path("results", "predictions"), 
  hive_style = TRUE, 
  unify_schemas = TRUE) |> 
  filter(
    modelid == run_config$model_id, 
    expected_type == expected_string, 
    year == 2021) 

## Get control IDs for matching

pred_id_control <- predY_ds |> 
  select(scenario, standid) |> 
  filter(scenario == "control") |> 
  collect() |> 
  pull(standid) |> 
  unique()

matches <- readRDS(file.path("data", "metso", "matched_pairs_utm.rds")) |> 
  mutate(in_pred_control = standid_matched_control %in% pred_id_control, row_id = row_number()) |> 
  filter(in_pred_control == 1) 

## Identify all IDs needed for the loop (N=100)

iterations <- 100
matches_subset <- matches[1:iterations, ]
needed_stands <- unique(c(matches_subset$standid_treated, matches_subset$standid_matched_control))
needed_utm    <- unique(c(matches_subset$UTM200_metso, matches_subset$UTM200_control))

## Batch collect 
predY <- predY_ds |> 
  filter(utm_zone %in% needed_utm) |> 
  filter(standid %in% needed_stands) |> 
  select(scenario, standid, posterior, all_of(name_spp)) |> 
  collect()


# 4. Run metrics

alpha <- 0.10
metrics <- list()

for(i in 1:iterations){
  
  # i <- 1
  # Filter
  matches_i <- matches_subset[i,]
  
  predY_i <- predY |> 
    filter(standid %in% c(matches_i$standid_treated, matches_i$standid_matched_control))
  
  # Matrix Wrangling
  
  predY_metso <- predY_i |> 
    filter(scenario == "metso") |> 
    arrange(posterior) |> 
    select(-scenario, -posterior, -standid) |> 
    as.matrix()
  
  predY_control <- predY_i |> 
    filter(scenario == "control") |> 
    arrange(posterior) |> 
    select(-scenario, -posterior, -standid) |> 
    as.matrix()
  
  # Vectorized Calculation
  
  if(nrow(predY_metso) == nrow(predY_control) && nrow(predY_metso) > 0) {
    metrics[[i]] <- calculate_metrics_vectorized(
      predY_metso_mat = predY_metso, 
      predY_control_mat = predY_control, 
      alpha = alpha
    )
  }
  
}
