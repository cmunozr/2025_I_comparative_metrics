# Experiments log

This file is the Narrative Journal. It is the logbook. The aim is to write down thoughts, hypotheses, plans, and conclusions in order to document the scientific reasoning behind of models that we have experimented. With this we can trace why we decided to create a specific model in the first place, what we learned from the model, and what your plan is for next running and models.


Structure of an entry

- Date
- Model ID at what is refering
- Miscelaneus notes
- Potential ideas to test (hypothesis)
- Action plan

## Entries

### September 19, 2025
#### Evaluating Model `fbs_M001`

This is a model constructed with all the species found in the Bakx list with data in the Finish bird Survey. Here, I selected the available variables from Shultz with a VIF anlysis. 

Covariates removed: average_stand_length, stand_basal_area, volume_total_wood,   tree_height_max.


**Miscelaneus notes:** The first full run of `fbs_M001` (see `Run_ID: fbs_M001_thin_100_samples_1000_chains_4`) showed a mix of good and poor MCMC convergence at 100 samples (in terms of ESS specifically).  Rho is really good. The PSRF mean values for beta, Gamma, and V are good although the maximum is to high (1.7951). Omega is mean PSRF > 1.1 in all random levels but sampleUnit. Alpha is also mean PSRF > 1.1 (1.528). Omega and Alpha has a maximum PSRF > 2.Actually, ALPHA AND RHO TRACE PLOTS LOOKS HORRIBLE.

**Ideas to test:** It is needed more thining on the sampling to improve the number of effective sample size.

**Action Plan:**
1.  Run `fbs_M001`with thin 1000 and samples 1000

### September 23, 2025
#### Evaluating Model `fbs_M001`thin 1000 and samples 1000

Model structure of `fbs_M001`

**Miscelaneus notes:** The second full run of `fbs_M001` (see `Run_ID: fbs_M001_thin_1000_samples_1000_chains_4`). Mean ESS is above 1000 in all paramaters. Rho keeps really good. Beta, gamma,  V, are good and the maximum PSRF < 1.1. Omega mean in all random levels is PSRF ~ 1.05 but omega maximum PSRF for vakio and sampleUnit is above 2 being the worst vakio. Alpha mean in vakio is 1.445 and maximum 1.997. Run the model with thin 1000 and samples 1000 improve the ESS and PSRF in almost all the parameters. ALPHA TRACE PLOTS LOOKS HORRIBLE.

**Ideas to test:** Following meetings with the team, I was advice to test the next settings: 

1) remove species under a threshold of summed abundance Tristan works with 35 https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.4559 ; 
2) latitud/longitud as variables in the model to reduce problems in the convergence of alpha and omega that use spatiall random levels, in this way we can modelling explicitly geography in the model; 
3) use altitude as covariate; 
4) remove occurrences above an altitud threshold; 
5) run a model with experst selection thorugh literature using the paper Wagenaar, L.F., Olsson, O., Stjernman, M. and Smith, H.G., 2025. A systematic meta-review: The relationship between forest structures and biodiversity in deciduous forests. Forest Ecology and Management, 596, p.123072; 
6) select variables using lasso; 
7) add climatic covariates
0) use shultz available variables without any filter-

**Action Plan:**
1. Set `fbs_M002`, that remove species under a threshold of summed abundance using 35
2. Run `fbs_M002`with thin 1000 and samples 1000

### September 30, 2025
#### Evaluating Model `fbs_M002`thin 1000 and samples 1000

Here, I remove species under a threshold of summed abundance Tristan works with 35 https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.4559. I selected covariates using VIF. To produce that effect S01 script was edited to create the `frequent_species` object in line 792 and then use to filter the occurrences of the species in subsections "occurr_sect_route".

Covariates removed: "average_stand_length"     "stand_basal_area"         "volume_total_wood"    "tree_height_max"    "canopy_cover_whole_stand"  and an alised factor "average_stand_length_gt250_mean" with "average_stand_length_gt200_mean"

**Miscelaneus notes:** Beta, Gamma, V are PSRF mean ~1.009 in ESS mean > 1000. Rho is 1.06. Omega kepps PSRF mean values around ~1.09. Here thw worst omega is sampleUnit with 3 but vakio takes 2.01 reducing it. Alpha PSRF is a bit better 1.385 in PSRF with a maximum of 2.264, making alpha worst than the another set. Potentially, in terms of alpha and omega random levels convergence, remove the rare species doesnt make much difference.

WARNING: Variable average_stand_diameter_gt15_mean was left out by accident. Potentially is needed to tun again in order to gain in systematic approach. However, removing this variable is not a hard problem. We continue the exploration of models.

**Action Plan:**
1. Set `fbs_M003`, it modells explicitly geography using longitud/latitud as covariate apart from the VIF selection of covariates, use an abundance threshold of 35.
2. Run `fbs_M003`with thin 1000 and samples 1000

### September 30, 2025
#### Creating Model `fbs_M003`thin 1000 and samples 1000

line 1248 add a chunck to calculate the centroid of each biotope and saved in `l2_xy` object, and line 1330 in S01 adds the coordinates of each biotope as covariate in XData using `cbind(l2_xy)`. Model running in Kingma.

### October 02, 2025
#### Creating Model `fbs_M004`thin 1000 and samples 1000

Elevation model was found in https://paituli.csc.fi/download.html but downloaded from the index https://www.nic.funet.fi/index/geodata/mml/dem10m/2019/. After extract and aggregation, now elevation is inside the pre-processed covariates. In that way elevation_mean is calculated automatically and added to the pool of covariates to use in XData. If you want to turn off elevation_mean as covariates just remove it in the line 1260 when `creating predictor_data` object of the same script. Model running in Orwell.









