#' DESCRIPTION:
#' Creates an Excel supplement with species-level survey data and model
#' performance for the selected presence--absence and abundance models.
#'
#' METHODOLOGY:
#' - Summarises cleaned observations from the selected fitted models.
#' - Extracts species-level fit, route-blocked CV and METSO hold-out metrics.
#' - Marks ratios as missing when their observed denominator is zero.
#' - Writes a README and three tables to a single Excel workbook.
#'
#' INPUTS:
#' - code/config_model.R: selected model runs.
#' - data/species_info_table.csv: species habitat classification.
#' - models/model_fit_master_all_metrics.csv: S05b model-fit results.
#' - models/{run}/fitted_{run}.rds: fitted models and response matrices.
#' - models/{run}/ho_metso/partition.rds: hold-out partition assignments.
#'
#' OUTPUT:
#' - results/supplementary/Supplementary_Tables.xlsx

library(dplyr)
library(openxlsx)

source(file.path("code", "config_model.R"))

model_path <- function(run, ...) file.path("models", run, ...)
read_model <- function(run) readRDS(model_path(run, paste0("fitted_", run, ".rds")))

pa_model <- read_model(selected_models$PA)
acp_model <- read_model(selected_models$aCp)
pa_holdout <- readRDS(model_path(selected_models$PA, "ho_metso", "partition.rds")) == 2
acp_holdout <- readRDS(model_path(selected_models$aCp, "ho_metso", "partition.rds")) == 2
species_info <- read.csv("data/species_info_table.csv", check.names = FALSE)
master_fit <- read.csv("models/model_fit_master_all_metrics.csv", check.names = FALSE)

species <- sort(colnames(pa_model$Y))
if (!setequal(species, colnames(acp_model$Y))) {
  stop("The selected models do not contain the same species.")
}

pa_y <- pa_model$Y[, species, drop = FALSE]
acp_y <- acp_model$Y[, species, drop = FALSE]
habitat <- species_info$Specialist[match(species, species_info$species)]

species_data <- data.frame(
  Species = species,
  `Habitat class` = habitat,
  `Recorded pairs` = colSums(acp_y, na.rm = TRUE),
  `Occupied route-years` = colSums(pa_y > 0, na.rm = TRUE),
  `Prevalence (%)` = 100 * colMeans(pa_y > 0, na.rm = TRUE),
  `Hold-out occupied route-years` = colSums(pa_y[pa_holdout, , drop = FALSE] > 0, na.rm = TRUE),
  check.names = FALSE
)

strategies <- c("Full_Model", "route_blocked_cv", "metso_holdout")
strategy_labels <- c("Explanatory", "Spatial CV", "METSO hold-out")

performance_table <- function(run, metrics) {
  fit <- master_fit[master_fit$Run_Name == run & master_fit$Strategy %in% strategies, ]
  index <- expand.grid(Strategy = strategies, Species = species, stringsAsFactors = FALSE)
  key <- paste(fit$Species, fit$Strategy, sep = "\r")
  wanted <- paste(index$Species, index$Strategy, sep = "\r")
  if (anyDuplicated(key) || anyNA(match(wanted, key))) {
    stop("Missing or duplicated species-level results for ", run, ". Run S05b first.")
  }
  fit <- fit[match(wanted, key), ]
  result <- data.frame(
    Species = index$Species,
    Evaluation = strategy_labels[match(index$Strategy, strategies)],
    fit[, unname(metrics), drop = FALSE],
    check.names = FALSE
  )
  names(result)[-(1:2)] <- names(metrics)
  result
}

pa_performance <- performance_table(selected_models$PA, c(
  "AUC" = "AUC", "Brier score" = "Brier", "Prevalence ratio" = "Prevalence_ratio",
  "Calibration slope" = "Calib_slope", "90% interval width" = "Mean_CI_width"
))
acp_performance <- performance_table(selected_models$aCp, c(
  "SR2" = "SR2", "MAE (pairs)" = "MAE", "Mean abundance ratio" = "Mean_ratio",
  "IQR dispersion ratio" = "IQR_ratio", "90% interval width (pairs)" = "Mean_CI_width"
))

holdout_zero <- colSums(pa_y[pa_holdout, , drop = FALSE] > 0, na.rm = TRUE) == 0
pa_performance$`Prevalence ratio`[
  pa_performance$Evaluation == "METSO hold-out" & holdout_zero[pa_performance$Species]
] <- NA_real_

has_observed_iqr <- function(y) {
  apply(y, 2, function(x) {
    x <- x[!is.na(x)]
    length(x) > 1 && IQR(x) > 0
  })
}
iqr_full <- has_observed_iqr(acp_y)
iqr_holdout <- has_observed_iqr(acp_y[acp_holdout, , drop = FALSE])
valid_iqr <- ifelse(
  acp_performance$Evaluation == "METSO hold-out",
  iqr_holdout[acp_performance$Species], iqr_full[acp_performance$Species]
)
acp_performance$`IQR dispersion ratio`[!valid_iqr] <- NA_real_

readme <- data.frame(
  Item = c("Purpose", "Table S1", "Table S2", "Table S3", "Species data",
           "Evaluations", "Hold-out", "Prevalence", "AUC", "Brier score",
           "Prevalence ratio", "Calibration slope", "SR2", "MAE", "Mean abundance ratio",
           "IQR dispersion ratio", "90% interval width", "Missing values", "Selected PA run", "Selected abundance run"),
  Description = c(
    "Species-level survey data and performance of the selected hurdle model components.",
    "Cleaned survey data and specialist/generalist habitat classification.",
    "Presence-absence performance by species and evaluation.",
    "Conditional-abundance performance by species and evaluation.",
    "Recorded pairs are summed after cleaning. Occupied route-years have at least one recorded pair. Prevalence is the percentage of all modelled FBS route-years occupied.",
    "Explanatory: full-data fit; Spatial CV: route-blocked cross-validation; METSO hold-out: targeted held-out routes.",
    paste0("Hold-out occupied route-years refer to the ", sum(pa_holdout), " route-years in the PA hold-out; total modelled route-years: ", nrow(pa_y), "."),
    "Prevalence is the proportion of modelled route-years with at least one recorded pair, expressed as a percentage.",
    "Area under the receiver operating characteristic curve; measures discrimination.",
    "Mean squared error of predicted occurrence probabilities; lower values are better.",
    "Predicted divided by observed prevalence; a value near one indicates agreement in average prevalence.",
    "Slope from observed occurrence regressed on the logit of predicted probability; one is the ideal value.",
    "Signed squared Spearman correlation between observed and predicted conditional abundance.",
    "Mean absolute error of conditional-abundance predictions, expressed in recorded pairs.",
    "Predicted mean divided by observed mean conditional abundance.",
    "Predicted interquartile range divided by observed interquartile range of conditional abundance.",
    "Mean width of the 90% posterior interval for the model's expected occurrence probability or conditional abundance, respectively.",
    "Blank means unavailable or undefined. Prevalence ratios with zero observed prevalence and IQR ratios with zero observed IQR are undefined.",
    selected_models$PA, selected_models$aCp
  ),
  check.names = FALSE
)

out_dir <- file.path("results", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workbook <- createWorkbook()
tables <- list("README" = readme, "Table S1 Species" = species_data,
               "Table S2 PA" = pa_performance, "Table S3 aCp" = acp_performance)
header_style <- createStyle(textDecoration = "bold", fgFill = "#DCEAF4", border = "bottom")
number_style <- createStyle(numFmt = "0.000")

for (name in names(tables)) {
  addWorksheet(workbook, name)
  writeData(workbook, name, tables[[name]], headerStyle = header_style, keepNA = FALSE)
  freezePane(workbook, name, firstRow = TRUE)
  if (name == "README") {
    setColWidths(workbook, name, cols = 1, widths = 25)
    setColWidths(workbook, name, cols = 2, widths = 110)
    addStyle(workbook, name, createStyle(wrapText = TRUE, valign = "top"),
             rows = 2:(nrow(readme) + 1), cols = 2, gridExpand = TRUE, stack = TRUE)
    setRowHeights(workbook, name, rows = 2:(nrow(readme) + 1), heights = 30)
  } else {
    addFilter(workbook, name, rows = 1, cols = seq_len(ncol(tables[[name]])))
    setColWidths(workbook, name, cols = 1, widths = 28)
    setColWidths(workbook, name, cols = 2:ncol(tables[[name]]), widths = "auto")
    if (name %in% c("Table S2 PA", "Table S3 aCp")) {
      addStyle(workbook, name, number_style, rows = 2:(nrow(tables[[name]]) + 1),
               cols = 3:ncol(tables[[name]]), gridExpand = TRUE, stack = TRUE)
    } else {
      addStyle(workbook, name, createStyle(numFmt = "0.0"), rows = 2:(nrow(species_data) + 1),
               cols = 5, gridExpand = TRUE, stack = TRUE)
    }
  }
}

output_path <- file.path(out_dir, "Supplementary_Tables.xlsx")
saveWorkbook(workbook, output_path, overwrite = TRUE)

# Remove unused drawing relationships created for worksheets without drawings.
workbook_dir <- tempfile("supplement_")
dir.create(workbook_dir)
unzip(output_path, exdir = workbook_dir)
relationship_dir <- file.path(workbook_dir, "xl", "worksheets", "_rels")
for (path in list.files(relationship_dir, pattern = "\\.rels$", full.names = TRUE)) {
  relationships <- xml2::read_xml(path)
  unused <- xml2::xml_find_all(relationships,
    "//*[local-name()='Relationship'][contains(@Type, '/drawing') or contains(@Type, '/vmlDrawing')]")
  xml2::xml_remove(unused)
  xml2::write_xml(relationships, path)
}
fixed_path <- tempfile(fileext = ".xlsx")
zip::zipr(fixed_path, files = list.files(workbook_dir, recursive = TRUE, all.files = TRUE), root = workbook_dir, mode = "mirror")
file.copy(fixed_path, output_path, overwrite = TRUE)
unlink(c(workbook_dir, fixed_path), recursive = TRUE)
message("Saved supplementary tables to: ", output_path)
