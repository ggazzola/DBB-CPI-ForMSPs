rm(list = ls())

setwd("/PathToYourWorkingDirectory") 
require(extendedForest)
source("VariableImportanceFunctions.R")


nRep = 100 
n = 100 

#ADD HERE THE PARAMETERS USED BY YOUR DATA-GENERATING FUNCTION
beta = # CHANGE ACCORDING TO THE IMPORTANCE THE VARIABLES IN YOUR SCENARIO HAVE ON THE FINAL OUTPUT
p = length(beta)


##############################################


TRUEScore = VINormalize(abs(beta)) 
TRUERank = rank(TRUEScore) 

standardize = T 

maxLev = 4
nTrees = 300 
corrTree = 0.15
mtr = ceiling(p/3)


impVPRFAllScore = NULL
impVPRFAllRank = NULL
perfIndivScore = perfIndivRank = NULL
datList= list()
set.seed(1) # YOU CAN USE THIS TO MAKE SURE EVERY TIME YOU REPEAT THE EXPERIMENT YOU GET THE SAME RESULTS
for (i in 1:nRep){

	dat = # CALL TO YOUR DATA GENERATING FUNCTION, WHICH SHOULD RETURN A MATRIX OR A DATAFRAME
	
	if (standardize)
		dat = as.data.frame(scale(dat)) 
	datList[[i]]=dat
}


for (i in 1:nRep) {
	dat = datList[[i]]	

	impVPRF = NULL

	RF <- randomForest(Y~., data = dat, maxLevel = maxLev, importance = TRUE, ntree = nTrees, corr.threshold = corrTree, mtry = mtr, keep.group=T)
	impVPRF = importance(RF)[,1]
	
	impVPRFAllScore = rbind(impVPRFAllScore, (impVPRF))
	impVPRFAllRank = rbind(impVPRFAllRank, rank(impVPRF))

	perfIndivScore[i] =  PerfImp(VINormalize(impVPRFAllScore[i,]), TRUEScore)
	perfIndivRank[i] = PerfImp(impVPRFAllRank[i,], TRUERank)

	print(paste(i/nRep*100, "% done"))

}

perfIndivScoreMean = 1/mean(perfIndivScore)
perfIndivRankMean= 1/mean(perfIndivRank)

VPRFScoreMean =(apply(impVPRFAllScore, 2, mean))
VPRFScoreSd =(apply(impVPRFAllScore, 2, sd))
VPRFRankMean =(apply(impVPRFAllRank, 2, mean))
VPRFRankSd =(apply(impVPRFAllRank, 2, sd))
perfMeanScore = 1/PerfImp(VPRFScoreMean, TRUEScore)
perfMeanRank = 1/PerfImp(VPRFRankMean, TRUERank)



res = list()	#THIS LIST CONTAINS ALL THE RESULTS WE NEED TO KEEP
res$Score = impVPRFAllScore
res$MeanScore = VPRFScoreMean
res$SdScore = VPRFScoreSd
res$Rank = impVPRFAllRank
res$MeanRank = VPRFRankMean
res$SdRank = VPRFRankSd

res$PerfIndivScore = 1/perfIndivScore
res$PerfIndivScoreMean= perfIndivScoreMean

res$PerfIndivRank = 1/perfIndivRank
res$PerfIndivRankMean= perfIndivRankMean

res$PerfMeanScore= perfMeanScore
res$PerfMeanRank= perfMeanRank

save("SomeFileName.RData", res) # SAVE THE RESULTS OF YOUR EXPERIMENTS FOR LATER ANALYSIS

