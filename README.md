
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
* [Python]() version 3.10 or upper
* This repository (2025_I_comparative_metrics)

### Libraries
R libraries required

```
matrixStats
Hmsc
tidyverse
data.table
terra
sf
maps
ape
taxize
remotes
phytools
openxlsx
tidyterra
viridis
lwgeom
ggrepel
gt
exactextractr) 
patchwork
dplyr
units
```

Python libraries required

```
tensorflow
hmsc-hpc
selenium
time
cdsapi
os
```

## How to run

We suggest running the routines step-by-step, following the order in the next flowchart. It relates 

IN CONSTRUCTION

## Authors and contact

## Main depeloving

* [Carlos Jair Muñoz Rodriguez](https://www.linkedin.com/in/carlos-jair-munoz/)
* [Ullrika Sahlin](https://www.cec.lu.se/ullrika-sahlin)

## Collaborators

* [Henrik Smith](https://www.cec.lu.se/henrik-smith)
* [Daniele Silvestro](https://bsse.ethz.ch/cevo/the-group/people/person-detail.dsilvestro.html)
* [Dario Shultz](https://www.ilr1.uni-bonn.de/en/about-us/staff/dario-schulz)
* [Aleksi Lehikoinen](https://researchportal.helsinki.fi/en/persons/aleksi-lehikoinen)


## Acknowledgment

* This development is supported by the [Mistra BIOPATH project](https://www.mistrabiopath.se/) - Aligning the Financial System with the Needs of Biodiversity. Also from the [Greenpole](https://greenpole.se/) an interdisciplinary research project generating evidence-based decision support for a sustainable use of Nordic forests. 

## License

This project is licensed under the MIT License - see the [License.md](LICENSE.md) file for details



