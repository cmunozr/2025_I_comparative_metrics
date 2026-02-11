
# Comparative Evaluation of Biodiversity Metrics for reporting Business Biodiversity Impact – an Assessment on Finnish Forest Birds 

This repository stores code and functions to automate the construction of Joint Species Distribution Models (J-SDM). These functions utilize Hmsc as the modelling engine. The Hmsc model provides predictions at the species-community level. This enables the spatial representation of biodiversity to analyse and evaluate biodiversity metrics’ performance and variation through environmental and spatial gradients while counting for uncertainties.

Current state: in development.

## Background

Biodiversity is declining rapidly, largely due to economic activities, posing significant risks to both the environment and financial sectors. Initiatives like BIOPATH are working to integrate biodiversity considerations into decision-making, which hinges on effective measurement. However, choosing the right biodiversity metrics is challenging, as there's a tension between business needs for simplicity and scientific needs for ecological complexity. The metric selected can significantly impact how biodiversity change is understood and addressed. This project aims to quantitatively compare various metrics to assess their accuracy in reflecting biodiversity states and trends under different forestry management scenarios, using birds as a focal taxon. We'll explore how well these metrics detect changes and how their sensitivity varies with sampling effort or spatial scale.

## Prerequisites

### Dependencies and files

Dependencies to install, choose the version depending on your operating system and version. For example, a windows 10 terminal with more than 4 gigabytes on memory RAM almost always has a 64 bit version of windows. Surf on the web in case of more information.

* [R](https://cran.r-project.org/mirrors.html) version 4.4.2 or upper
* [RStudio](https://www.rstudio.com/products/rstudio/download/#download) version 1.4 or upper
* [Python](https://www.python.org/downloads/) version 3.10 or upper (required for GPU acceleration)
* This repository (2025_I_comparative_metrics)

### Libraries

#### R Libraries

Standard CRAN packages used across scripts:

| Package | Version |
|---|---|
| matrixStats | 1.5.0 |
| Hmsc | 3.3.7 |
| tidyverse | 2.0.0 |
| data.table | 1.17.8 |
| terra | 1.8.86 |
| sf | 1.0.23 |
| maps | 3.4.3 |
| ape | 5.8.1 |
| taxize | 0.10.0 |
| remotes | 2.5.0 |
| phytools | 2.5.2 |
| openxlsx | 4.2.8.1 |
| tidyterra | 0.7.2 |
| viridis | 0.6.5 |
| lwgeom | 0.2.14 |
| ggrepel | 0.9.6 |
| gt | 1.1.0 |
| exactextractr | 0.10.0 |
| patchwork | 1.3.2 |
| dplyr | 1.1.4 |
| units | 1.0.0 |
| car | 3.1.3 |
| arrow | 23.0.0 |
| zip | 2.3.3 |
| geometry | 0.5.2 |
| jsonify | 1.2.3 |
| here | 1.0.2 |
| FD | 1.0.12.3 |

GitHub/Custom packages (installation required via devtools or remotes):

```
# install_github("eliotmiller/clootl")
clootl
# install_github("luomus/finbif")
finbif
```

#### Python Libraries

These are required for the hmsc-hpc engine (see code/gpu_instructions.txt for setup):

```
tensorflow[and-cuda]
tensorflow-probability
hmsc-hpc
time
cdsapi
os
```

## Documentation

For a deep dive into the code architecture, function parameters, and detailed methodology, please visit the full documentation:
**[Project Wiki & Documentation](https://deepwiki.com/cmunozr/2025_I_comparative_metrics)**

This README provides a quick-start guide. For detailed function references (e.g., how `mean_species_abundance` is calculated), refer to the Wiki.

## How to run

The workflow is organized into sequential scripts (prefixed S01, S02, etc.) located in the `code/` directory.

### 1. Data Preparation & Model Definition
* **S01_define_models_finland.qmd**:
  * Loads Finnish Bird Survey (FBS) data, harmonizes taxonomy (FinBIF, ITIS), and filters by biotope.
  * Prepares covariates (Forest inventory, climate) and traits/phylogeny.
  * Constructs the unfitted Hmsc model and saves it to `models/unfitted_[ID].RData`.

### 2. Model Fitting (GPU Accelerated)
Fitting is performed using Python/TensorFlow for speed.

* **S02a_fit_models_gpu.R**: Reads the unfitted model and generates the necessary Python scripts (`.sh` files) and JSON configurations for the HPC/GPU environment.
* **Execution**: Run the generated shell scripts on your GPU server (refer to `code/gpu_instructions.txt` for specific commands like `nohup ... &`).
* **S02b_fit_models_gpu_import_posterior.R**: Once MCMC sampling is complete, this script imports the posterior samples back into the R Hmsc object.

### 3. Diagnostics & Evaluation
* **S03 Series**: Visualizes MCMC convergence (trace plots, effective sample size).
* **S04 Series**: Calculates model fit indices (WAIC, AUC, Tjur's R2) using partition (CV) or internal evaluation.
* **S05 Series**: Generates plots summarizing model explanatory power.

### 4. Ecological Inference
* **S06_show_parameter_estimates.R**: Visualizes Beta parameters (environmental responses) and Gamma parameters (traits).
* **S07_visualize_responses.R**: Plots gradient responses for specific species or community traits.

### 5. Scenario Analysis & Metric Comparison
* **S08_predict_scenarios_utmSplit.R**: Generates posterior predictions for specific forestry scenarios (e.g., METSO conservation vs. Business-As-Usual).
* **S09_run_diversity_metrics.R**:
  * Calculates taxonomic and functional diversity metrics on the predicted communities.
  * Uses parallel processing (batching) to handle large posterior arrays.
* **S10 Series**: Analyzes and visualizes the distribution of these metrics across scenarios.

### Data Requirements
Before running the modeling scripts (S01+), you must ensure the covariate and metso data is available. 
* **Automatic Download:** You can use the provide download scripts (see next) to get the required datasets (Metsakeskus, Luke, etc.).
* **Note:** Ensure you have sufficient disk space as defined in the source files.

## Utilities and Project Architecture

In addition to the main workflow scripts (`S01`–`S10`), this repository contains `_utilities` scripts. These serve two distinct purposes:

### 1. Function Libraries (Backend)
These files contain core mathematical functions, model definitions, and prediction logic.
* **Do not run these directly.** They are automatically `source()`'d by the main analysis scripts (e.g., `S02`, `S08`, `S09`) as needed.

* `_utilities_diversity_metrics_functions.R` (Calculates MSA, Rao's Q, etc.)
* `_utilities_hmsc_gpu.R` (HPC/GPU configurations)
* `_utilities_transform_covariates.R`

### 2. Data Preparation Tools (Run-Once)
These are standalone scripts used to fetch raw data or construct specific datasets (e.g., the "Metso" conservation areas). Run these **only** if you are building the project data from scratch.

#### Phase A: Initial Setup (Run before S01)
* **Covariate Acquisition:**
    * `_utilities_download_covariates.R`
    * `_utilities_download_covariates_climatic.py`
    * `_utilities_download_metsakeskus.R`
* **Covariate Pre-processing:**
    * `_utilities_preprocess_climatic_covariates.R`
* **Study Design & Treatment Units (Metso):**
    * `_utilities_metso_trt_control_setup.R`
    * `_utilities_metso_merge_metsakeskus.R`
    * `_utilities_metso_trt_control_labelbydist.R`
* **Extracts values for Metso sites**
    * `_utilities_fetch_covariates.R` 

#### Phase B: Scenario Construction (Run before S08)
Before running the scenario predictions in `S08`, you must generate the predictive covariate matrices (`XData`) for the different management scenarios.

* `_utilities_metso_calc_XData.R`

#### Phase C: Diagnostics & Benchmarking (Optional)
These scripts are useful for testing the prediction engine before committing to the full prediction in `S08`.

* `_utilities_predict_scenarios_benchmark.R`
  * *Purpose:* Benchmark processing time using different batch sizes to optimize performance.
* `_utilities_predict_scenarios_replicates.R`
  * *Purpose:* Ensure the replication and consistency of the prediction generation process.


## Authors and contact

## Main developers

* [Carlos Jair Muñoz Rodriguez](https://www.linkedin.com/in/carlos-jair-munoz/)
* [Ullrika Sahlin](https://www.cec.lu.se/ullrika-sahlin)

## Collaborators

* [Henrik Smith](https://www.cec.lu.se/henrik-smith)
* [Daniele Silvestro](https://bsse.ethz.ch/cevo/the-group/people/person-detail.dsilvestro.html)
* [Dario Shultz](https://www.ilr1.uni-bonn.de/en/about-us/staff/dario-schulz)
* [Aleksi Lehikoinen](https://researchportal.helsinki.fi/en/persons/aleksi-lehikoinen)

## Acknowledgment

This work has been carried out with support from the [Mistra BIOPATH project](https://www.mistrabiopath.se/) — Aligning the Financial System with the Needs of Biodiversity. Additional support comes from [Greenpole](https://greenpole.se/), an interdisciplinary research effort generating evidence‑based decision tools for the sustainable management of Nordic forests.

## License

This project is licensed under the MIT License - see the [License.md](LICENSE.md) file for details
