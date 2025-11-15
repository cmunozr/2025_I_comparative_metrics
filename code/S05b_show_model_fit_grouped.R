# --- 1. Load Libraries and Configuration ---
library(tidyverse)
library(dplyr)
library(stringr)
library(colorspace) 

# --- 2. Setup Directories and Initialize Data Lists ---
model_setup <- list.dirs("models", recursive = FALSE) |>
  basename() |>
  stringr::str_subset("^fbs_M")

labels <- c("cv", "ho")

convergence_folder <- file.path("models", model_setup, labels)
convergence_folder[sapply(convergence_folder, file.exists)]

output_dir <- file.path("models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message(paste0("--- Starting comparison for all models runned ---"))

# Initialize lists to store aggregated data from all runs
psrf.beta.all <- list(); ess.beta.all <- list()
psrf.gamma.all <- list(); ess.gamma.all <- list()
psrf.rho.all <- list(); ess.rho.all <- list()
psrf.V.all <- list(); ess.V.all <- list()
psrf.alpha.all <- list(); ess.alpha.all <- list() 
psrf.omega.vakio.all <- list(); ess.omega.vakio.all <- list()
psrf.omega.year.all <- list(); ess.omega.year.all <- list()
psrf.omega.sampleUnit.all <- list(); ess.omega.sampleUnit.all <- list()

# statistical function to summarize
sum_func <- median
sum_name <- "Median"