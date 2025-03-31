# THIS SCRIPT MAKES PREDICTIONS OVER OVER SPACE FOR HMSC MODELS FOR THE BIRD EXAMPLE (SECTION 11.1) OF THE BOOK
# Ovaskainen, O. and Abrego, N. 2020. Joint Species Distribution Modelling - With Applications in R. Cambridge University Press.

# Set the base directory

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Fitted model and prediction grid.

#	OUTPUT. Spatial predictions illustrated as various plots.
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


library(Hmsc)

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nst = length(thin_list)
nChains = 4

# We first read in the model object from the longest RUN
for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){break}
}
load(filename)

m = models[[1]]

# Here we perform spatial predictions. To do so, we proceed as in the single-species case sudy in Chapter 4.
# Thus, we import a grid of spatial coordinates, habitat types and climatic conditions for 1000 locations.
# We then apply the `prepareGradient` function to these data to prepare a spatial gradient for which the
# predictions are to be made
# While in the book the predictions are presented for 10000 prediction locations, we recommend running this
# script for 1000 prediction locations to make the running time faster. To choose which one to do,
# read either the file "grid_1000.csv" to "grid_10000.csv"

grid = read.csv(file.path(dataDir, "grid_1000.csv"), stringsAsFactors=TRUE)

# Let's look at what is included in the grid

head(grid)

# We see that there  are the same environmental variables as in the data used to fit the model,
# as well as the spatial coordinates. This is what we need to make the predictions.

# The habitat type "Ma" is present in the prediction data but it is not part of training data,
# so we can't make predictions for it

grid = droplevels(subset(grid,!(Habitat=="Ma")))

# We next construct the objects xy.grid and XData.grid that have the coordinates and
# environmental predictors, named similarly (hab and clim) as for the original data matrix (see m$XData) 

xy.grid = as.matrix(cbind(grid$x,grid$y))
XData.grid = data.frame(hab=grid$Habitat, clim=grid$AprMay, stringsAsFactors = TRUE)

# We next use the prepareGradient function to convert the environmental and spatial
# predictors into a format that can be used as input for the predict function

Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(Route=xy.grid))
#
# We are now ready to compute the posterior predictive distribution
# takes couple of minues to compute it for 1000 locations

if(TRUE){
  nParallel=2
  predY = predict(m, Gradient=Gradient, expected = TRUE, nParallel=nParallel)
  
  # Note that we used expected = TRUE to predict occurrence probabilities (e.g. 0.2) instead of occurrences (0 or 1) 
  # Note also that if you have very large prediction grid, you can use the predictEtaMean = TRUE option
  # to speed up the computations
  # predY = predict(m, Gradient=Gradient, predictEtaMean = TRUE, expected = TRUE)
  
  # Let's explore the prediction object.
  
  class(predY)
  
  # It is a list... 
  
  length(predY)
  
  # ...of length 1000, if you fitted four chains with 250 samples from each
  
  dim(predY[[1]])
  
  # Each prediction is a matrix with dimensions 951 x 50, as there are 951 prediction locations and 50 species.
  
  head(predY[[1]])
  
  # Each matrix is filled in with occurrence probabilities
  # We may simply by ignoring parameter uncertainty and just looking at 
  # the posterior mean prediction. 
  
  EpredY=Reduce("+",predY)/length(predY)
  dim(EpredY)
  head(EpredY)
  
  # EpredY is a 951 x 50 matrix of posterior mean occurrence probabilities
  
  # Let us save this object so we do not need to recompute it if continuing with this work
  save(EpredY, file = "results/predictions.RData")
}
load("results/predictions.RData")

# The next step is to post-process the predictions to those community features
# that we wish to illustrate over the prediction space. With the script below,
# we derive from the predictions the occurrence probability of C. monedula (species number 50),
# the species richness, and community-weighted mean traits.
# We also include data on habitat type and climatic conditions to the dataframe mapData that
# includes all the information we need to visualize the predictions as maps

Cm = EpredY[,50]
S=rowSums(EpredY)
CWM = (EpredY%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
xy = grid[,1:2]
H = XData.grid$hab
C = grid$AprMay
RCP10 = kmeans(EpredY, 10)
RCP10.cluster = as.factor(RCP10$cluster)
mapData=data.frame(xy,C,S,Cm,CWM,H,RCP10.cluster,stringsAsFactors=TRUE)

# We will use the ggplot function from the ggplot2 package, so let's load the data
library(ggplot2)

# We first plot variation in the habitat and climatic conditions on which the predictions are based on.

ggplot(data = mapData, aes(x=x, y=y, color=H))+geom_point(size=1) + ggtitle("Habitat") + scale_color_discrete() + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=C))+geom_point(size=1) + ggtitle("Climate") + scale_color_gradient(low="blue", high="red") + coord_equal()

# We then exemplify prediction for one focal species, here C. monedula over Finland

ggplot(data = mapData, aes(x=x, y=y, color=Cm))+geom_point(size=1) + ggtitle(expression(italic("Corvus monedula")))+ scale_color_gradient(low="blue", high="red") + coord_equal()

# This prediction is reassuringly very similar to that based on the single-species model of Chapter 5.7
# We next plot predicted species richness, which is highest in Southern Finland

ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=1) + ggtitle("Species richness")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the proportion of resident species, also highest in Southern Finland 

ggplot(data = mapData, aes(x=x, y=y, color=MigrationR))+geom_point(size=1) + ggtitle("Proportion of resident species")+ scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot the community-weighted mean log-transformed body size, which is highest in Northern Finland

ggplot(data = mapData, aes(x=x, y=y, color=LogMass))+geom_point(size=1) + ggtitle("Mean log-transformed body mass") + scale_color_gradient(low="blue", high="red") + coord_equal()

# We next plot regions of common profile
ggplot(data = mapData, aes(x=x, y=y, color=RCP10.cluster))+geom_point(size=1) + ggtitle("RCP10") + scale_color_discrete() + coord_equal()
