
nParallel = 4 #Default: nParallel = nChains
#nParallel = 1 #Set to 1 to disable parallel computing

localDir = getwd()
modelDir = file.path(localDir, "models")

library(Hmsc)

load(file = file.path(modelDir,"unfitted_models.RData"))
nm <- length(models)

if(.Platform$OS.type == "windows"){
  # LAPTOP
  samples_list = c(5,250,250)
  thin_list <-  c(1,1,10)
}else{
  samples_list = c(5, 250,250, 250 )
  thin_list <-  c(1, 100,1000, 10000)
}

nChains <-  4
if(is.null(nParallel)) nParallel <- nChains
Lst = 1

library(doParallel)
registerDoParallel(cores = 4)

while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = ""))
  if(file.exists(filename)){
    print("model had been fitted already")
  } else {
    print(date())
    for (mi in 1:nm) {
      print(paste0("model = ",names(models)[mi]))
      m = models[[mi]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel) 
      models[[mi]] = m
    }
    save(models,file=filename)
  }
  Lst = Lst + 1
}
