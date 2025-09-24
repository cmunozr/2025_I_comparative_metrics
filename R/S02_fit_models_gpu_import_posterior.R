library(jsonify)
library(Hmsc)

chains <- 4
thin <- 1000
samples <- 1000
transient <- 500000
dir.models <- "models"

# DEFAULT RUN
# 
# load("models/unfitted_models.RData")
# m <- models$swedishBirds_forestry_abundance_small_ngpp
# gpu_run <- from_json(readRDS("models/swedishBirds_forestry_abundance_GPU_burn_5000_thin_100_samples_1000_chains_4_post_file.rds_temp")[[1]])
# postList <-  gpu_run[1:4]
# 
# fitTF <-  importPosteriorFromHPC(m = m, postList = postList, nSamples = 250, thin = 1, transient = 125, alignPost = T)

#------------------

## RUN BY CHAINS

load("models/unfitted_fbsF_001.RData")
nm_to_save <- "fitted_fbsF_002.RData"

nm <- names(models)

for(i in 1:length(nm)){
  # i <- 1
  nm_i <- nm[i]
  
  root <- paste0(nm_i, "_GPU", "_thin_", thin, "_samples_", samples, "_chains_", chains)
  
  chain_python <- paste0("0", 0:(chains-1))
  
  branch <- paste0(root, "_post_chain", chain_python, "_file.rds")
  
  complete_path <- file.path(dir.models, root, branch)
  
  runned_logic <- lapply(X = complete_path, FUN = function(X){file.exists(X)}) |> unlist()
  
  if(sum(runned_logic) < chains){
    warning(paste0("MODEL: ", nm_i, ". RUNNED CHAINS ARE LESS THAN THE NUMBER OF EXPECTED CHAINS (", chains, "). POST FILE WILL BE CONSTRUCTED WITH: ", sum(runned_logic), " CHAINS"))
  }
  
  runned_path <- complete_path[runned_logic]
  
  runned_list <- list() 
  
  for(r in 1:length(runned_path)){
    #r <- 1
    runned_r <- runned_path[r]
    runned_list[[r]] <- from_json(readRDS(runned_r)[[1]])[[1]]
  }
  
  runned_post <-  importPosteriorFromHPC(m = models[[i]], postList = runned_list, nSamples = samples, thin = thin, transient = transient, alignPost = T)
  
  models[[i]] <- runned_post
  
}

save(models, file = file.path(dir.models, nm_to_save))

rm(list=ls());gc()





