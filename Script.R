setwd("/Users/gianluca/Desktop/AllDesktop/Stuff/Rutcor/Research/Samsung/2014/VariableImportance/Code")
require(relaimpo)
require(pls)
require(CORElearn)
require(extendedForest)
source("VariableImportanceFunctions.R")
source("RelWeights.R")
source("VIP.R")

p =5
n = 100
rho = 0
beta = c(6,0,3,1,3)
noiseSd = 10#22
standardize = T



datGen = DataGenGromping(n, p, rho, beta, noiseSd)
dat = datGen$dat
if (standardize)
	dat = as.data.frame(scale(dat))
# Dominance Analysis - LMG
linMod = lm(Y~., data = dat)
impGDW = calc.relimp(linMod, type = "lmg", rela=T)
impGDW = impGDW@lmg

# Relative Weights
impRW = relweights(linMod, col = "lightgrey")
impRW = t(impRW)
row.names(impRW)="impRW"
#VIP
pls <- plsr(Y~., p, data = dat, validation = "none", method = "oscorespls") #10 is the number of components? %SCALE FIRST?
impVIP<- VIP(pls)[p,]

# Variable Permutation with RF, %maxlevel 0 gives you marginal permutation, K considers the first K splits (that is 2^K groups)
# if mtry = 1, then splitting variable chosen completely at random (the splitting criterion has no effect), maybe a good choice otherwise bias toward variables that are important for the Y
RF <- randomForest(Y~., data = dat, maxLevel = 4, importance = TRUE, ntree = 500, corr.threshold = 0.5, mtry = 1)
impVPRF = importance(RF)[,1]

# RELIEF
impREL <- attrEval(Y~., data = dat, estimator="RReliefFexpRank", kNearestExpRank = n/2, ReliefIterations=0) #kNearestEqual, #kNearestExpRank

# simple regression

impREG = abs(linMod$coefficients[-1])

# expected score
impTRUE = abs(beta)
results =rbind(impGDW, impRW, impVIP, impVPRF, impREL, impREG, impTRUE)

(t(apply(results, 1, rank)))

round(t(apply(results, 1, VINormalize)),2)
