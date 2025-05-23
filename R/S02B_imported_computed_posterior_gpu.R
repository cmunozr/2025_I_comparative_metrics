library(jsonify)
library(Hmsc)

load("models/unfitted_models.RData")
m <- models$swedishBirds_forestry_abundance_small_ngpp
gpu_run <- from_json(readRDS("models/swedishBirds_forestry_abundance_GPU_burn_5000_thin_100_samples_1000_chains_4_post_file.rds_temp")[[1]])
postList <-  gpu_run[1:4]

fitTF <-  importPosteriorFromHPC(m = m, postList = postList, nSamples = 250, thin = 1, transient = 125, alignPost = T)

