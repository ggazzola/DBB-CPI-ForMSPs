
DataGenSlides = function(n, standardize) { # old data generating model, you can disregard this
	alpha = sqrt(1-0.3^2)
	epsilon =NULL;
	outX = outY =NULL
	for (j in 1:n) {
		for (i in 1:5) 
			epsilon[i] = rnorm(1)
		
		x1 = 10 + epsilon[1];
		x2 = 10 + 0.3*epsilon[1] + alpha*epsilon[2]
		x3 = 10 + 0.3*epsilon[1] + 0.5604*alpha*epsilon[2] + 0.8282*alpha*epsilon[3]
		x4 = -8 + x1 - 0.5*x2 + 0.3*x3 + 0.5*epsilon[4]
		x5 = -5 + 0.5*x1 + x2 + 0.5*epsilon[5]
		x6 = 6*rnorm(1)
		x7 = 7*rnorm(1)
		x8 = 8*rnorm(1)
		x9 = 9*rnorm(1)
		x10 = 10*rnorm(1)
		# x0 = 1
		outX = rbind(outX, c(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))
		y = - 8 + x1 + 0.5*x2 + 0.3*x3 + 0.5*epsilon[1]
		outY = c(outY, y)
	}
	
	if (standardize) {
		outX = scale(outX)
		outY = scale(outY)
	}	
	out = data.frame(X = outX, Y = t(t(outY)))
	return(out)
}

DataGenGromping = function(n, p, rho, beta, noiseSd, extraInputs) { #data generating model as described in my report PrelResSimulNov29.doc
	covMat = CovMatGen(p, rho, rhoIdxDist = T)
	datX = mvrnorm(n , rep(0,p), covMat)
	if (length(extraInputs)>0) {
		for (i in extraInputs)
			datX[,i]= i*rnorm(n)
	}	
	datY = datX %*%beta + rnorm(n)*noiseSd
	res = DatRes(datX, datY)
}

DataGenSquare = function(n, p , rho, rhoIdxDist, betaLin, betaSq, noiseSd, extraInputs) {
	# same correlation scheme as DataGenGromping, but input variables are squared
	covMat = CovMatGen(p, rho, rhoIdxDist)
	datX = mvrnorm(n , rep(0,p), covMat)
	if (length(extraInputs)>0) {
		for (i in extraInputs)
			datX[,i]= runif(n)
	}
	datY = datX%*%betaLin + datX^2 %*%betaSq + rnorm(n)*noiseSd
	res = DatRes(datX, datY)
}


DataGenBlock = function(n, blockP, blockNum, meanVect, rho, rhoIdxDist, noiseSd, linear) { # blocked input variable data generating model
	#n : number of observations
	#blockP: number of input variables in each block
	# blockNum: number of blocks
	# meanVect: mean vector of all input variables, across all blocks
	# rho: controls variable correlation within blocks
	# rhoIdxDist: set = T if want covariance structure as in DataGenGromp within block; set = F if want all vaiables within block to have the same pairwise correlation
	# noiseSd: variance of output noise
	# linear: set = T if linear generating model wanted, set = F if exponential generating model wanted
	covMat = matrix(0, blockP*blockNum, blockP*blockNum)
	rhoSeq = seq(0, rho, length.out = blockNum) #this is used to make the within block correlations increase linearly from the first to the last block 
												# (e.g. see rho_b in my "ReportDec1.doc" file
	cnt = 1
	for (i in 1:blockNum){
		currRho = rhoSeq[i]
		currBlockCov = CovMatGen(blockP, currRho, rhoIdxDist)
		covMat[cnt:(cnt+blockP-1), cnt:(cnt+blockP-1)] = currBlockCov
		cnt = cnt + blockP
	}
	
	datX = mvrnorm(n, meanVect, covMat)
	if (linear) {
		datY = datX %*%beta + rnorm(n)*noiseSd
	} else {
		datY = exp(datX %*%beta)/(1+exp(datX %*%beta))+ rnorm(n)*noiseSd
	}	
	res = DatRes(datX, datY)
}

CovMatGen = function(p, rho, rhoIdxDist) { # creation of covariance matrix. Set rhoIdxDist input to TRUE to reproduce the scheme explained in PrelResSimulNov29.doc for Sigma
	covMat = matrix(0, p, p) # assuming input variables have variance 1
	for (i in 1:p) {
		for (j in 1:p) {
			if (i==j) {
				covMat[i,j] = 1 
			} else {
				if (rhoIdxDist)
					covMat[i,j] = rho^(abs(i-j))
				else
					covMat[i,j] = rho
			}
		}
	}
	return(covMat)
}


VINormalize = function(VIVect) { # normalization of variable importance scores so that they sum to 100
	if (min(VIVect)<0)
		VIVect[VIVect<0]=0
		#VIVect=VIVect+abs(min(VIVect)) ## maybe
	if (sum(VIVect)>0) {
		multiplier = 100/sum(VIVect)
	} else{
		multiplier = 0
	}
	VIVect = multiplier * VIVect
	return(VIVect)
}

#VINormalize = function(VIVect) { # normalization of variable importance scores so that they sum to 100
#	VIVect = (VIVect-min(VIVect))/(max(VIVect)-min(VIVect))
#		return(VIVect)
#}

DatRes = function(datX, datY) { # estimation of linear regression model and its R^2 from a data set with inputs datX and output datY
	dat = data.frame(X = datX, Y = datY)
	linFit = lm(Y~., data=dat)
	rSquared = summary(linFit)$r.squared
	res = list(dat = dat, rSquared = rSquared)
}


PerfImp = function(obsImpVect, trueImpVect, normType = "L2") { #used for the computation of P1 and P2
	diffVect = obsImpVect - trueImpVect
	if (normType == "L2") 
		nrm = sum(diffVect^2)
	else
		nrm = sum(abs(diffVect))
	
	res = nrm
	#res = 1/nrm 
	return(res)
}