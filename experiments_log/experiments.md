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


**Miscelaneus notes:** This model showed a mix of good and poor MCMC convergence at 100 samples (in terms of ESS specifically).  Rho is really good. The PSRF mean values for beta, Gamma, and V are good although the maximum is to high (1.7951). Omega is mean PSRF > 1.1 in all random levels but sampleUnit. Alpha is also mean PSRF > 1.1 (1.528). Omega and Alpha has a maximum PSRF > 2.Actually, ALPHA AND RHO TRACE PLOTS LOOKS HORRIBLE.

**Ideas to test:** It is needed more thining on the sampling to improve the number of effective sample size.

**Action Plan:**
1.  Run `fbs_M001`with thin 1000 and samples 1000

### September 23, 2025
#### Evaluating Model `fbs_M001` thin 1000 and samples 1000

Model structure of `fbs_M001`

**Miscelaneus notes:** Mean ESS is above 1000 in all paramaters. Rho keeps really good. Beta, gamma,  V, are good and the maximum PSRF < 1.1. Omega mean in all random levels is PSRF ~ 1.05 but omega maximum PSRF for vakio and sampleUnit is above 2 being the worst vakio. Alpha mean in vakio is 1.445 and maximum 1.997. Run the model with thin 1000 and samples 1000 improve the ESS and PSRF in almost all the parameters. ALPHA TRACE PLOTS LOOKS HORRIBLE.

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
#### Evaluating Model `fbs_M002` thin 1000 and samples 1000

Here, I remove species under a threshold of summed abundance Tristan works with 35 https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.4559. I selected covariates using VIF. To produce that effect S01 script was edited to create the `frequent_species` object in line 792 and then use to filter the occurrences of the species in subsections "occurr_sect_route".

Covariates removed: "average_stand_length"     "stand_basal_area"         "volume_total_wood"    "tree_height_max"    "canopy_cover_whole_stand"  and an alised factor "average_stand_length_gt250_mean" with "average_stand_length_gt200_mean"

**Miscelaneus notes:** Beta, Gamma, V are PSRF mean ~1.009 in ESS mean > 1000. Rho is 1.06. Omega kepps PSRF mean values around ~1.09. Here thw worst omega is sampleUnit with 3 but vakio takes 2.01 reducing it. Alpha PSRF is a bit better 1.385 in PSRF with a maximum of 2.264, making alpha worst than the another set. Potentially, in terms of alpha and omega random levels convergence, remove the rare species doesnt make much difference.

WARNING: Variable average_stand_diameter_gt15_mean was left out by accident. Potentially is needed to tun again in order to gain in systematic approach. However, removing this variable is not a hard problem. We continue the exploration of models.

**Action Plan:**
1. Set `fbs_M003`, it modells explicitly geography using longitud/latitud as covariate apart from the VIF selection of covariates, use an abundance threshold of 35.
2. Run `fbs_M003`with thin 1000 and samples 1000

### September 30, 2025
#### Creating Model `fbs_M003` thin 1000 and samples 1000

line 1248 add a chunck to calculate the centroid of each biotope and saved in `l2_xy` object, and line 1330 in S01 adds the coordinates of each biotope as covariate in XData using `cbind(l2_xy)`. Model running in Kingma.

### October 02, 2025
#### Creating Model `fbs_M004` thin 1000 and samples 1000

Elevation model was found in https://paituli.csc.fi/download.html but downloaded from the index https://www.nic.funet.fi/index/geodata/mml/dem10m/2019/. After extract and aggregation, now elevation is inside the pre-processed covariates. In that way elevation_mean is calculated automatically and added to the pool of covariates to use in XData. If you want to turn off elevation_mean as covariates just remove it in the line 1260 when `creating predictor_data` object of the same script. Model running in Orwell.

### October 06, 2025
#### Checking runnings of `fbs_M003` and `fbs_M004`

Issue in output folder in the server side. Need to run everything again. Script S02 was refactored. 

### October 08, 2025
#### Checking runnings of `fbs_M003` and `fbs_M004`

- `fbs_M003` memory error for chains 00 and 01. Chains 02 and 03 completed. When the server was accesed the run for chain 03 was launched before finish chain 01. As the models are running in automatic, we need to check the bash file.
- `fbs_M004` running as expected. Chain 00 and 01 ready. Chain 02 and 03 sampling.

### October 09, 2025
#### Evaluating Model `fbs_M004` thin 1000 and samples 1000

As we have problems running model `fbs_M003`, we are commenting first this model, the `fbs_M004` model. 

This model works using the threshold of summed abundance based in Baxk work. Also, use coordinates (centroid of each biotope polygon) as a covariate and elevation mean. Selected covariates using VIF. See entrie on October 02 for more info.

**Miscelaneus notes:** Parameters Beta, Gamma, and V keep a mean PSRF near to 1 and mean ESS of above 1000. Rho has a bad value, 2.8. When checking the actual values of Rho the chains show high values (0-1). However, each one moved in different values like chain 1 = 0.97, chain 2 = 0.98, chain 3 = 0.99, chain 4 = 1. *In this way, althoug the PSRF has a poor measure, can we say that actually the model has converged?* (note trace plots are not good neither "stuck"). 

"PSRF compares the variance within each chain to the variance between chains. If the chains have converged, these variances should be similar, and the PSRF value will approach 1"

About parameter Omega at all three parameters according to the mean is good: 1.02-1.07. PSRF max of Omega varies between 2.3-4.0. Some min ESS are exceptionally low as for random level year (27.73) or sampleUnit (105.9).

About the alpha parameter convergence and ESS are quite good.

*Another question that arises is: what is the meaning, in the context of niche modeling, of using coordinates and altitude as covariates? In the tropics, altitude is closely related to temperature, but this relationship may differ in temperate or boreal regions.*

### October 10, 2025
#### Evaluating Model `fbs_M003` thin 1000 and samples 1000

This model works using the threshold of summed abundance based in Baxk work. Also, use coordinates (centroid of each biotope polygon) as a covariate. Selected covariates using VIF. See entrie on September 30 for more info.

**Miscelaneus notes:** Parameters Beta, Gamma, and V keep a mean PSRF near to 1 and mean ESS of above 1000. Rho has good value, 1.01. All chains are 1 for Rho. *Is there a convergence? The samples arent really moving, or maybe PSRF is calculated through all the values of the chain?* I dont think

About parameter Omega at all three parameters according to the mean is good: 1.01-1.03. PSRF max of Omega varies between 1.5-3.2. Some min ESS are exceptionally low as for random level year (36.94) or sampleUnit (89.42).

About the alpha parameter convergence (mean psrf = 1.04, max = 1.3) and ESS are good.

*Same question as the model `fbs_M004`: what is the meaning, in the context of niche modeling, of using coordinates as covariates?*

**Action Plan:**
1. Using this simpler model, comparing with `fbs_M004`, can we reach a better convergence in the maximum values if I take more samples? Run `fbs_M003`with thin 1000 and samples 2000.

### October 10, 2025
#### Creating Model `fbs_M003` thin 1000 and samples 2000

Model `fbs_M003` reach a correct convergence and distribution as shown with PSRF and ESS of almost all parameters. However, I think is problematic that the maximum value is above the threshold of a good model. Lets try to take more samples and see if the proceses reach a better convergence and the minimum ESS is better for the omega parameters and V.

### October 16, 2025
#### Evaluating Model `fbs_M003` thin 1000 and samples 2000

**Miscelaneus notes:** This model has a good convergence and effective sample size in Beta, Gamma and V. However, Rho convergence is worst with 3.86 and actually values on the omega values are worst for temporal and unit level around 1.1. The alpha has good convergence. Minimun values for omega random levels are really low as last attempts as the maximum values are equally big. I dont think that having a longer run could improve the convergence of the model, actually in general terms is a bit worst, maybe because the sampler is visiting several Rho has the same behaviour as in the model `fbs_M004` thin 1000 and samples 1000.

It is not worthy to spend 6 days fitting a model with 2000 samples. Also, the potential failing of chains is higher as the writing of the rds using pyreadr issue is not fixed yet, although there is a potential solution. https://github.com/hmsc-r/hmsc-hpc/pull/2

### October 16, 2025
#### General 

**Miscelaneus notes:** in entrie mark with September 30, 2025 I dismissed the accidental deletion of the variable: average_stand_diameter_gt15_mean. However, as I continue the comparisson between model `fbs_M002` and the rest is being hard to complete as the basic structure is not the same as the others. Although is needed to say that when we remove species the structure of the correlation between variables 'vary' and in this way the models are not completely comparable.

Aditionally, when re-visiting model `fbs_M001` thin 1000 samples 1000 its discrepancy with the model `fbs_M003` thin 1000 samples 1000 is in the Alpha parameter, the former with mean 1.445 and the latter with 1.071. What if we can make converge a model with the complete set of species using the spatial fixed effects. Although, here comes again the need for a justification of use coordinates as a covariate.

**Action Plan:**

1. Run model `fbs_M005` which will be structurally similar to the foundational model `fbs_M001` but using coordinates as a covariate.
2. Re run model `fbs_M002` with name `fbs_M002.1` making sure that the basic variables are the same or left out for correlation issues and not of a carelessness handle from the researcher.

### October 16, 2025
#### Creating Model `fbs_M002.1` thin 1000 and samples 1000

S02 script for deefine models was transdormed to deactivate the binding of object l2_xy deactivated and avoiding coordinates a fixed effect. Check entry from September 30, 2025. Equally, deleted elevation_mean from XData as it is calculated automatically from the utilities fetch from covariates script. line ~1334 look for the object XData.



