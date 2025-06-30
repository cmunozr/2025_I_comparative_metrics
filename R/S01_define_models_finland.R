library(openxlsx)
library(dplyr)

#------------

stand_template <- read.xlsx("slu/ImportTemplate_ver3.xlsx", sheet = "StandData")
treatment_template <- read.xlsx("slu/ImportTemplate_ver3.xlsx", sheet = "TreatmentProposals")

#------------

coord <- read.csv("slu/input_heureka_coordSampleUnits.csv")
forest <- read.csv("slu/input_heureka_forest.csv")


