
##	Efficient number of species

##	Beta: turnover


if (!requireNamespace("betapart", quietly = TRUE)) install.packages("betapart")
library(betapart)

a <- betapart.core(logic_predY[1:100,,1]*1)

beta.pair(x = logic_predY[,,1]*1, index.family = "sorensen")



multiple_dissimilarity(pred.object = logic_predY[,,1:2], aggregation.factor = 2, richness = a){
  
  
  
  
}

##	Beta: nestedness 

ex <- vegan::beta.div(x = predY[1:100, , 1])
plot(ex)