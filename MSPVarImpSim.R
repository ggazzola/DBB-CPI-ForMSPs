require(relaimpo)	#package 'relaimpo' should be installed prior to runnning this code
require(gtools)
source("CYPaths.R") #the file "CYPaths.R" should be in the R working directory
set.seed(1)
# MSP DATA GENERATION PARAMETERS

n = 1000	#number of observations
Q = 8	#total number of stages
minX = 1	#range of X variables is [minX, maxX]
maxX = 10
sigmaEpsVect = rep(1,Q) #standard deviation of error term for the Q stages
betaZeroVect = rep(1,Q) #intercept of true model for the Q stages
depType = "full"	#type of inter-stage dependency (can be "sequential or "full")
standardize = T	# standardize data along each X and Y variable?

# WHICH CONTRIBUTION SHOULD BE COMPUTED?

s = 1	#stage to compute the contribution from
t = Q	#stage to compute the contribution to

XMat = YMat = EpsMat = NULL
depList = impList = cList =list()


# GENERATE DATA

for (i in 1:Q) {
	XMat = cbind(XMat, runif(n, minX, maxX)) #Sampling values for X variables in the i-th stage
	EpsMat = cbind(EpsMat, rnorm(n, 0, sigmaEpsVect[i]))	#Sampling values for error variables in the i-th stage
	
	depList[[i]] = list(dep=NULL, Xcoef = NULL, Ycoef = NULL)	#Adjacency matrix: entry (i,j)=1 if stage i depends on stage j; 
																#(i,j)=0 otherwise. 1/0 values are assigned just below
	if (standardize)
		XMat[,i] = as.matrix(scale(XMat[,i]))
		
	if(i>1) {
		if (depType=="sequential") {
			depList[[i]]$dep = i-1 # sequential dependencies
		} else {	
			depList[[i]]$dep = 1:(i-1) # full dependencies
		}
	}
	
	tmp = rep(0,Q) 
	if (i>1)
		tmp[depList[[i]]$dep] = rep(1,length(depList[[i]]$dep))# seq(1, 1/length(depList[[i]]$dep), length.out=length(depList[[i]]$dep)) # a sequence of decreasing coefficients
	
	depList[[i]]$Ycoef = tmp	#coefficients for the Y regressors in the i-th true model (as set just above in 'tmp')
	depList[[i]]$Xcoef = 1#0.5	#coefficient for the X regressor in the i-th true model
	
	XComp = as.matrix(XMat[,i])%*%depList[[i]]$Xcoef #multiplication of X regressors by their respective coefficients
	if (i>1) {
		nonZeroYCoef = depList[[i]]$Ycoef
		nonZeroYCoef = nonZeroYCoef[nonZeroYCoef!=0]
		if(standardize){
			YComp = as.matrix(scale(YMat[,depList[[i]]$dep]))%*% nonZeroYCoef #inner product of Y regressors by their respective coefficients
		} else{
			YComp = as.matrix(YMat[,depList[[i]]$dep])%*% nonZeroYCoef #inner product of Y regressors by their respective coefficients
			
		}	
	} else{
		YComp = 0	# since there are no Y regressors in the first stage, there is no inner product to compute 
	}
	YMat = cbind(YMat, XComp + YComp + EpsMat[,i] + betaZeroVect[i]) #computing output of i-th stage for all observations by summing all inputs
}

colnames(XMat) = paste("X",1:Q, sep="")
colnames(YMat) = paste("Y",1:Q, sep="")
colnames(EpsMat) = paste("Eps",1:Q, sep="")


# DISPLAY CORRELATIONS IN ARTIFICIAL MSP DATA

allVarMat = as.data.frame(cbind(XMat, EpsMat, YMat))
round(cor(allVarMat),2)


# ASSESS VARIABLE IMPORTANCE WITHIN EACH INDIVIDUAL STAGE

for (i in 1:Q){
	impList[[i]] = list(formula = NULL, model = NULL, rsq = NULL, relImp = NULL)
	relImp  = rep(0,Q+1)

	if(i >1 ) {
		YInput = depList[[i]]$dep
		form = as.formula(paste("Y", i, "~", "X", i, "+", paste("Y", YInput, sep="",collapse="+"), sep="")) # definition of linear model formula for the i-th stage
		impList[[i]]$model = lm(form, data = allVarMat)	# estimation of linear model for the i-th stage
		tmp = calc.relimp(impList[[i]]$model, type = "lmg", rela=T)@lmg #computation of regressors' relative importance using general dominance weights ("lmg" method in package relaimpo)
		relImp[1] = tmp[1] ######## assuming X is always present (we could do the same with Y, in case we don't know dependencies b/w Y's a priori)
		relImp[YInput + 1] = tmp[-1]
		impList[[i]]$relImp = relImp	#the first element of 'relImp' is the relative importance of the X regressor; the remaining elements are the relative importance of the Y regressors
										#note that if a given Y variable does not appear as regressor in the linear model for the i-th stage, its corresponding element in relImp is = 0
										#note also that the sum of the elements of relImp is equal to 1
	} else{
		form = as.formula(paste("Y", i, "~", "X", i, sep=""))	# same as above for first stage
		relImp[1] =  1
		impList[[i]]$model = lm(form, data = allVarMat)
		impList[[i]]$relImp = relImp
	}	
	names(impList[[i]]$relImp) = c(paste("X", i, sep=""), paste("Y", 1:Q, sep=""))
	impList[[i]]$formula = form
	impList[[i]]$rsq = summary(impList[[i]]$model)$r.squared # coefficient of determination of linear model for the i-th stage

}

# BREAK DOWN RELATIVE IMPORTANCE BY STAGE-INHERENT INPUT VARIABLE (X), STAGE-INHERENT NOISE VARIABLE (EPSILON), 
# AND ANY INPUT VARIABLES DEFINED BY PREVIOUS STAGES' OUTPUTS (Y'S);
# CREATE WEIGHTED MSP ADJACENCY MATRIX 

cYMatFromTo = NULL
for (i in 1:Q){
	currImp = impList[[i]]
	numInput = length(depList[[i]]$dep)
	cList[[i]] = list(cX=NULL, cEps=NULL, cYVect=NULL)  #the elements of 'cList[[i]]', as computed below, correspond to the share of the i-th stage output variance
														#explained by X regressor ($cX), noise term ($cEps), and Y regressors ($cYVect)
	 											
	cList[[i]]$cX = currImp$rsq*currImp$relImp[1]	
	cList[[i]]$cEps = 1 - currImp$rsq
	cList[[i]]$cYVect = currImp$rsq*currImp$relImp[-1]		
	cYMatFromTo = cbind(cYMatFromTo, cList[[i]]$cYVect)	#this is a weighted version of the adjacency matrix computed above, such that entry (i,j)  
	 													#corresponds to the share of variance of Yi explained by Yj
}
colnames(cYMatFromTo) = row.names(cYMatFromTo)

cS  = as.numeric(cList[[s]]$cEps + cList[[s]]$cX) #computation of c_s (as per formula (14), page 13, March 25 report) for stage 's' 
result = CYPaths(s,t,cYMatFromTo, cS) #computation of c_{s,t} (as per formula (14), page 13, March 25 report), using function 'CYPaths' 
totalContribution = result$totContrib
contributionByPath = result$pathList
# FIND ALL PATHS FROM CONTRIBUTING STAGE TO CONTRIBUTED STAGE

CYPaths = function(start, end, weightedAdjMat, cStart){
	#On the network induced by weighted adjacency matrix 'weightedAdjMat', this function finds all paths from stage 'start' to stage 'end'
	#For each such path, it computes the share of variance of Y{end} (output of stage 'end') explained by Y{start} (output of stage 'start') along the path,
	#corresponding to the product π shown in formula (14), page 13, March 25 report
	#'cStart' corresponds to c_s in the same formula
	
	cProdPathList = cProdPathListRes = list()
	path = c(start, end)
	
	if (start < end) {
		cProdPath = cStart * weightedAdjMat[start, end]
	} else {
		cProdPath = cStart 
	}
	
	cProdPathList[[1]]=list(length = 1, paths = path, cProd = cProdPath)	
	
	if (start<end) {
		cDirect = weightedAdjMat[start, end]
		tmp = start:end
		tmpMid = tmp[c(-1, -length(tmp))]
		tmpMidN  = length(tmpMid)
		if (tmpMidN > 0){
			for (i in 1:tmpMidN){
				currMidNodes = combinations(tmpMidN, i, tmpMid)
				cProdPathVect = pathMat = NULL
				cnt = 1
				for (j in 1:nrow(currMidNodes)){
					path = c(start, currMidNodes[j,], end)
					cProdPath = cStart
					for (k in 1:(length(path)-1)){
						cProdPath = cProdPath*weightedAdjMat[path[k], path[k+1]]
						if(cProdPath==0)
							break
					}
					pathMat = rbind(pathMat, path)
					rownames(pathMat)= NULL
					cProdPathVect = c(cProdPathVect, cProdPath )
					cnt = cnt + 1
				}
				cProdPathList[[i+1]] = list(length = length(path)-1, paths = pathMat, cProd = cProdPathVect)
			}
		}
	}
	
	sumContrib = 0
	for (i in 1:length(cProdPathList)) {
		sumContrib = sumContrib + sum(cProdPathList[[i]]$cProd)
	}

	res =list(pathList = cProdPathList, totContrib = sumContrib) 	#pathList[[i]] corresponds to paths of length i from stage 'start' to stage 'end' ($length); 
	return(res)														# The sequence of nodes defining the possible paths of length i are shown as rows in $paths
}																	# The relative contribution of 'start' to 'end' along the j-th path in $paths is given by the 
																	# j-th element in $cProd. 
																	# If the j-th element of $cProd is equal to 0, then the j-th path listed in $paths does not 
																	# exist in the MSP network (or the contribution along that path is exactly 0).

																	#totContrib gives the overall contribution of 'start' to 'end' along all paths from 'start' to 'end'
													


