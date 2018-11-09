rm(list = ls())

setwd("/Users/gianluca/Desktop/AllDesktop/Stuff/Rutcor/Research/Samsung/2014/VariableImportance/Code") 
require(relaimpo)	
require(pls)
require(CORElearn)
source("RelWeights.R")
require(MASS)
source("VariableImportanceFunctions.R")
source("VIP.R")
source("ExtractSplits.R")

allMethods = F
GGG = T

if(GGG) {
	require(extendedForestGGGFinal) 
	pow = 7 #### PARAMETER
} else {
	require(extendedForest)
}

#INSTRUCTIONS:
#REPORT RESULTS IN ONE EXCEL FILE PER EXPERIMENT (FOR WRITING DATA TO A FILE, YOU CAN USE, E.G., R FUNCTIONS write.table or write.csv )
#IN SHEET ONE, INCLUDE THE EXPERIMENTAL PARAMETERS ANY COMMENTS/OBSERVATIONS, IF YOU FIND OUT SOMETHING INTERESTING
#IN SHEET TWO INCLUDE ALL XXXSCOREMEAN, XXXSCORESD, XXXRANKMEAN, XXXRANKSD
#IN SHEET THREE INCLUDE ALL perfMeanScore, perfMeanRank, perfIndivScore, perfIndivRank, AND HIGHLIGHT THE BEST RESULT (HIGHEST VALUE WITHIN EACH perfXXXX)


#FOR THE TIME BEING, VARY JUST ONE PARAMETER AT A TIME, KEEPING ALL OTHERS FIXED 

#UN-COMMENT IF USING SCENARIO 1 
#REPORT ALL OF THESE EXPERIMENTAL PARAMETERS#########
nRep = 100 # DO NOT DECREASE, UNLESS RUNNING ONE EXPERIMENT WITH THE CURRENT EXPERIMENTAL PARAMETERS TAKES MORE THAN A FEW MINUTES
extraInputs = 6:8 #VARY, MAKE SURE IT'S CONSISTENT WITH BETA VECTOR ABOVE: THE REDUNDANT VARIABLES INDEXED BY THIS VECTOR SHOULD BE ASSOCIATED WITH ZERO BETA COEFFICIENTS
beta = c(c(5,5,5,5,5), rep(0, length(extraInputs))) #VARY (CHANGE (1) COEFFICIENT VALUES, E.G., DIFFERENT MAGNITUDES, WITH/WITHOUT TIES BETWEEN VARIABLES, WITH/WITHOUT ZERO COEFFICIENTS, ETC.; (2) NUMBER OF INPUT VARIALES, I.E., SHORTER OR LONGER BETA VECTOR)
n = 100 # VARY (E.G., 20, 1000)
rho = 0.75 #VARY
noiseSd = 3 #VARY
p = length(beta) 
##############################################


#UN-COMMENT IF USING SCENARIO 2 
#REPORT ALL OF THESE EXPERIMENTAL PARAMETERS#######
#nRep = 100 # DO NOT DECREASE, UNLESS RUNNING ONE EXPERIMENT WITH THE CURRENT EXPERIMENTAL PARAMETERS TAKES MORE THAN A FEW MINUTES
#blockP = 3 #VARY
#blockNum = 3 #VARY
#whichBlock = 3 #VARY (WHICH BLOCK CONTAINS RELEVANT INPUT VARIABLES)
#rhoIdxDist = FALSE #VARY (FALSE GIVES INTRA-BLOCK CORRELATION MATRICES WITH CONSTANT PAIRWISE CORRELATIONS, TRUE GIVES INTRA-BLOCK CORRELATION MATRICES WITH PAIRWISE CORRELATIONS THAT VARY DEPENDING ON THE VARIABLE INDEX, AS IN SCENARIO 1)
#n = 100 # VARY (E.G., 20, 1000)
#betaOne = 1 #VARY
#beta = rep(0, blockP*blockNum)				#VARY THIS AND, ACCORDINGLY, THE LINE BELOW (DEFINITION OF BETA VECTOR). THE CURRENT SETTING ASSIGNS A COEFFICIENT = betaOne TO THE FIRST VARIABLE OF THE whichBlock-TH BLOCK. ALL OTHER COEFFICIENTS ARE ZERO
#beta[(whichBlock-1)*blockP+1] = betaOne
#noiseSd = 2 #VARY
#rho = 0.7 #VARY
#meanVect = rep(0, blockNum*blockP)
#p = length(beta) 
#linear = TRUE 
##############################################


TRUEScore = VINormalize(abs(beta)) # the true input variable score, given the generating model, normalized so that scores sum to 100
TRUERank = rank(TRUEScore) # the true input variable rank, given the generating model

standardize = T # standardize the data 

#Parameters for randomForest function below. Check its help documentation for details, but you can probably leave them as they are
maxLev = 5
nTrees = 300 #the higher this number, the more stable the results of variable permutation with random forest should be
corrTree = 0.15
mtr = ceiling(p/3)


impGDWAllScore = impRWAllScore  = impVIPAllScore  = impVPRFAllScore  = impRELAllScore  = impREGAllScore  = NULL
impGDWAllRank = impRWAllRank = impVIPAllRank = impVPRFAllRank = impRELAllRank = impREGAllRank = NULL
perfIndivScore = perfIndivRank = matrix(0, nRep,6)
perfMeanScore = perfMeanRank= matrix(0, 1,6)
datList= list()
set.seed(1) ###########
for (i in 1:nRep){
	datGen = DataGenGromping(n, p, rho, beta, noiseSd, extraInputs) #UN-COMMENT IF USING SCENARIO 1 
	#datGen = DataGenBlock(n, blockP, blockNum, meanVect, rho, rhoIdxDist, noiseSd, linear) #UN-COMMENT IF USING SCENARIO 2
	dat = datGen$dat
	if (standardize)
		dat = as.data.frame(scale(dat)) #data standardization
	
	datList[[i]]=dat
}

rm(.Random.seed)

for (i in 1:nRep) {
	dat = datList[[i]]	
	
	if(allMethods) {
		# General dominance weights method
		linMod = lm(Y~., data = dat)
		impGDW = calc.relimp(linMod, type = "lmg", rela=T)
		impGDW = impGDW@lmg
		impGDWAllScore = rbind(impGDWAllScore, (impGDW)) #
		impGDWAllRank = rbind(impGDWAllRank, rank(impGDW)) #


		# Relative weights method
		impRW = relweights(linMod, col = "lightgrey")
		impRW = t(impRW)
		row.names(impRW)="impRW"
		impRWAllScore = rbind(impRWAllScore, (impRW))
		impRWAllRank = rbind(impGDWAllRank, rank(impRW))
	
		# Variable importance in projection method
		pls <- plsr(Y~., p, data = dat, validation = "none", method = "oscorespls")  
		impVIP<- VIP(pls)[p,]
		impVIPAllScore = rbind(impVIPAllScore, (impVIP))
		impVIPAllRank = rbind(impVIPAllRank, rank(impVIP))
	}
	# Variable permutation with random forests method
	impVPRF = NULL

	if (GGG) {
		corDat = abs(cor(dat))
		corDat = corDat[1:p,1:p] ### assuming Y is p+1-st column of dat
		freqSplits = list()
		for (j in 1:p) {
			probVect = (corDat[j,])^pow
			probVect[j] =mean(probVect[-j]) 
			probVect = probVect/sum(probVect) #Assigning variable of interest a probability that = the average of the others
			RF <- randomForest(Y~., data = dat, maxLevel = maxLev, importance = TRUE, ntree = nTrees,  corr.threshold = corrTree, mtry = mtr,  ### keep corrThreshold>0 in Final?
			keep.group=T, prob=probVect) #varToScore=j,
			freqSplits[[j]] = apply(SplitFreq(RF,maxLev,p), 2, mean)

			impVPRF[j] = importance(RF)[j,1] #if importance is zero, it means that that variable was never used to split, altough it was always included in the candidate variables
		}										#if variable has low correlation with all others, then it's gonna be only predictor:. what kind of model and results for that case, esp. if variable is irrelevant?
	} else{
		RF <- randomForest(Y~., data = dat, maxLevel = maxLev, importance = TRUE, ntree = nTrees, corr.threshold = corrTree, mtry = mtr, keep.group=T)
		impVPRF = importance(RF)[,1]
		freqSplits = apply(SplitFreq(RF,maxLev,p), 2, mean)

	}
	impVPRFAllScore = rbind(impVPRFAllScore, (impVPRF))
	impVPRFAllRank = rbind(impVPRFAllRank, rank(impVPRF))

	if(allMethods){
		# R-Relief method
		impREL <- attrEval(Y~., data = dat, estimator="RReliefFexpRank", kNearestExpRank = ceiling(n/2), ReliefIterations=0)
		impRELAllScore = rbind(impRELAllScore, (impREL))
		impRELAllRank = rbind(impRELAllRank, rank(impREL))
	
		# Absolute regression coefficient method
		impREG =linMod$coefficients[-1]
		impREGAllScore = rbind(impREGAllScore, (abs(impREG)))
		impREGAllRank = rbind(impREGAllRank, rank(abs(impREG)))

		# here we build the averages that define performance measure P2 
		perfIndivScore[i,1] =  PerfImp(VINormalize(impGDWAllScore[i,]), TRUEScore)	
		perfIndivScore[i,2] = PerfImp(VINormalize(impRWAllScore[i,]), TRUEScore)
		perfIndivScore[i,3] = PerfImp(VINormalize(impVIPAllScore[i,]), TRUEScore)
	}	
	perfIndivScore[i,4] =  PerfImp(VINormalize(impVPRFAllScore[i,]), TRUEScore)
	
	if(allMethods){
	
		perfIndivScore[i,5] = PerfImp(VINormalize(impRELAllScore[i,]), TRUEScore)
		perfIndivScore[i, 6] = PerfImp(VINormalize(impREGAllScore[i,]), TRUEScore)

		perfIndivRank[i,1] = PerfImp(impGDWAllRank[i,], TRUERank)
		perfIndivRank[i,2]	= PerfImp(impRWAllRank[i,], TRUERank)
		perfIndivRank[i,3] = PerfImp(impVIPAllRank[i,], TRUERank)
	}	
	
	perfIndivRank[i,4] = PerfImp(impVPRFAllRank[i,], TRUERank)
	
	if(allMethods){
		perfIndivRank[i,5] = PerfImp(impRELAllRank[i,], TRUERank)
		perfIndivRank[i,6] = PerfImp(impREGAllRank[i,], TRUERank)
	}
	print(paste(i/nRep*100, "% done"))

}

perfIndivScoreMean = 1/apply(perfIndivScore,2,mean)
perfIndivRankMean= 1/apply(perfIndivRank,2,mean)

### REPORT THESE##########################
VPRFScoreMean =(apply(impVPRFAllScore, 2, mean))

if(allMethods){
	GWDScoreMean =(apply(impGDWAllScore, 2, mean)) #normalization ?
	RWScoreMean =(apply(impRWAllScore, 2, mean))	
	VIPScoreMean =(apply(impVIPAllScore, 2, mean))
	RELScoreMean =(apply(impRELAllScore, 2, mean))
	REGScoreMean =(apply(impREGAllScore, 2, mean))
##############################################

### REPORT THESE##########################
	GWDScoreSd =(apply(impGDWAllScore, 2, sd))
	RWScoreSd =(apply(impRWAllScore, 2, sd))	
	VIPScoreSd =(apply(impVIPAllScore, 2, sd))
}	
VPRFScoreSd =(apply(impVPRFAllScore, 2, sd))

if(allMethods){

	RELScoreSd =(apply(impRELAllScore, 2, sd))
	REGScoreSd =(apply(impREGAllScore, 2, sd))
##############################################

### REPORT THESE##########################
	GWDRankMean =(apply(impGDWAllRank, 2, mean))
	RWRankMean =(apply(impRWAllRank, 2, mean))	
	VIPRankMean =(apply(impVIPAllRank, 2, mean))
}
	
VPRFRankMean =(apply(impVPRFAllRank, 2, mean))

if(allMethods){

	RELRankMean =(apply(impRELAllRank, 2, mean))
	REGRankMean =(apply(impREGAllRank, 2, mean))
##############################################

### REPORT THESE##########################
	GWDRankSd =(apply(impGDWAllRank, 2, sd))
	RWRankSd =(apply(impRWAllRank, 2, sd))	
	VIPRankSd =(apply(impVIPAllRank, 2, sd))
}	
VPRFRankSd =(apply(impVPRFAllRank, 2, sd))

if(allMethods) {
	RELRankSd =(apply(impRELAllRank, 2, sd))
	REGRankSd =(apply(impREGAllRank, 2, sd))
##############################################

	perfMeanScore[1] = 1/PerfImp(GWDScoreMean, TRUEScore)
	perfMeanScore[2] = 1/PerfImp(RWScoreMean, TRUEScore)
	perfMeanScore[3] = 1/PerfImp(VIPScoreMean, TRUEScore)
}
perfMeanScore[4] = 1/PerfImp(VPRFScoreMean, TRUEScore)
if(allMethods){
	perfMeanScore[5] = 1/PerfImp(RELScoreMean, TRUEScore)
	perfMeanScore[6] = 1/PerfImp(REGScoreMean, TRUEScore)

	perfMeanRank[1] = 1/PerfImp(GWDRankMean, TRUERank)
	perfMeanRank[2] = 1/PerfImp(RWRankMean, TRUERank)
	perfMeanRank[3] = 1/PerfImp(VIPRankMean, TRUERank)
}	
perfMeanRank[4] = 1/PerfImp(VPRFRankMean, TRUERank)

if(allMethods){

	perfMeanRank[5] = 1/PerfImp(RELRankMean, TRUERank)
	perfMeanRank[6] = 1/PerfImp(REGRankMean, TRUERank)

colnames(perfIndivScore)=colnames(perfMeanScore)=colnames(perfIndivRank)=colnames(perfMeanRank)=c("GDW", "RW", "VIP", "VPRF", "REL", "REG")
}
### REPORT THESE##########################
#MULTIPLY perfMeanScore BY WHATEVER POWER OF 10 IT TAKES TO MAKE ALL ITS ENTRIES  > 1 
#E.G, IF perfMeanScore = (0.09163895 0.09766494 0.07308208 0.1052411 0.1072595 1.95515), MULTIPLY IT BY 10^2 TO OBTAIN 
# (9.163895 9.766494 7.308208 10.52411 10.72595 195.515). REPORT THIS LAST VECTOR AND THE POWER OF 10 YOU USED
#SAME THING FOR perfMeanRank
#perfMeanScore = as.data.frame(perfMeanScore)	
#perfMeanRank = as.data.frame(perfMeanRank)	
##############################################

res = list()
res$Score = impVPRFAllScore
res$MeanScore = VPRFScoreMean
res$SdScore = VPRFScoreSd
res$Rank = impVPRFAllRank

res$PerfMeanScore= perfMeanScore
res$PerfIndivScoreMean= perfIndivScoreMean
res$MeanRank = VPRFRankMean
res$SdRank = VPRFRankSd
res$PerfMeanRank= perfMeanRank
res$PerfIndivRankMean= perfIndivRankMean
if(GGG) save.image("ResC.RData") else save.image("ResS.RData")

#apply(SplitFreq(RF,maxLev,p), 2, mean)
