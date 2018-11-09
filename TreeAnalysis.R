
CondPermImpGGG <- function(varOfInterest, condVarVect, maxDepth, dat, rfPredObj, rfCondObj, catIdxVect, catValueList, nonDirectional = F){
	# computes conditional variable importance with one RF for modeling/prediction 
	# 	and a separate one for data partitioning
	
	# varOfInterest: an integer, corresponding to the index of the variable to compute the importance of 
	# condVarVect: a vector of integers, corresponding to the indices of the variables to condition upon
	# maxDepth: an integer, corresponding to the depth at which every tree (and corresponding rules) 
	#	of the partitioning RF is to be pruned
	# dat: a matrix/data frame of X-Y observations
	# rfPredObj: a RF object to be used for modeling/prediction (trained on dat)
	# rfPredObj: a RF object to be used for partitioning (trained on dat)
	
	totVarNum = ncol(dat)-1
	numTree = rfPredObj$ntree
	mseIncreaseVect = rep(0, numTree)
	for(i in 1:numTree) {
		treeIdx = i
		condTree = as.data.frame(getTree(rfCondObj,treeIdx))
		getCondsOut = GetConds(condTree, catIdxVect, catValueList)

		oobLogical = !as.logical(rfPredObj$inbag[,treeIdx]) 
		datOob = dat[oobLogical,] 

		prunePartOut = PrunePart(getCondsOut, depth = maxDepth, nonDirectional)
		if(nonDirectional) {
			extractPartOut = ExtractPart(prunePartOut, totVarNum)
			extractCombinedRulesOut = ExtractCombinedRules(extractPartOut, condVarVect)
		} else {
			extractCombinedRulesOut = prunePartOut
		}
		if (is.null(extractCombinedRulesOut)) { 
			partitionDataOut =list()
			partitionDataOut[[1]] = list(datPart = datOob, n = nrow(datOob))
		} else {
			partitionDataOut = PartitionData(extractCombinedRulesOut, datOob) 
		}
		removeEmptyPartitionsOut = RemoveEmptyPartitions(partitionDataOut)
		cc = 0
		for(kk in 1:length(removeEmptyPartitionsOut)) {
			cc = cc+removeEmptyPartitionsOut[[kk]]$n
		}
		if(cc!=sum(oobLogical))
			stop(paste("Tree", i))
		
		if(length(removeEmptyPartitionsOut)>0){ # if at least one partition has some data
			permuteVarOfInterestOut  = PermuteVarOfInterest(removeEmptyPartitionsOut, varOfInterest)
			rawPermPredictionsOut = RawPermPredictions(permuteVarOfInterestOut, rfPredObj, treeIdx) 
			mseIncreaseOut = MSEIncrease(rawPermPredictionsOut)
			mseIncreaseVect[i] = mseIncreaseOut
		}
	}
	res = mean(mseIncreaseVect)
	return(res)
}


MSEIncrease <- function(rawPermPredictionsOut) {
	#returns the MSE increase after permutation of the variable of interest
	#computed in all data subsets as per rawPermPredictionsOut
	numPart = length(rawPermPredictionsOut)
	permErrIncrease = 0
	res = rep(0, numPart)
	for(i in 1:numPart){
		rawErr = mean((rawPermPredictionsOut[[i]]$trueY - rawPermPredictionsOut[[i]]$rawPred)^2)
		permErr = mean((rawPermPredictionsOut[[i]]$trueY - rawPermPredictionsOut[[i]]$permPred)^2)
		permErrIncrease = permErrIncrease + (permErr-rawErr)
	}
	return(permErrIncrease)
}

RawPermPredictions <- function(permuteVarOfInterestOut, rfPredObj, treeIdx) {
	# for each raw data  subset and its permutation along the variable of interest
	# as contained in permuteVarOfInterestOut, it uses the treeIdx-th tree of rfObj 
	# to predict the output of the raw and permuted data set
	# rfPredObj should be a randomForest object predictor
	# treeIdx should be an integer
	
	res = list()
	for(i in 1:length(permuteVarOfInterestOut)) {
		datPartX = permuteVarOfInterestOut[[i]]$datPartX
		datPartXPerm = permuteVarOfInterestOut[[i]]$datPartXPerm
		datPartY = permuteVarOfInterestOut[[i]]$datPartY
		rawPred = predict(rfPredObj, VectToMat(datPartX), predict.all=T)$individual[,treeIdx]
		permPred = predict(rfPredObj, VectToMat(datPartXPerm), predict.all=T)$individual[,treeIdx]
		res[[i]]=list(rawPred = rawPred, permPred = permPred, trueY = datPartY)
	}
	return(res)
}

PermuteVarOfInterest <- function(removeEmptyPartitionsOut, varOfInterest) {
	# permutes variable of interest within each data subset given by removeEmptyPartitionsOut
	# varOfInterest should be an integer from 1 to p
	res = list()
	for (i in 1:length(removeEmptyPartitionsOut)) {
		currDat = removeEmptyPartitionsOut[[i]]$datPart
		currN  = removeEmptyPartitionsOut[[i]]$n
		colNames = colnames(currDat)
		whichY = which(colNames == "Y")
		currX = VectToMat(currDat[,-whichY])
		currY = currDat[,whichY]

		permVarOfInterest  = currX[,varOfInterest] 
		permVarOfInterest = permVarOfInterest[sample(currN)]

		permutedCurrDat = currX
		permutedCurrDat[, varOfInterest] = permVarOfInterest
		
		res[[i]] = list(datPartX = currX, datPartXPerm = permutedCurrDat, datPartY = currY)
	}
	return(res)
}

RemoveEmptyPartitions <- function(partitionDataOut){
	# remove any element of the partitionDataOut list
	# corresponding to partitions with no oob data
	res = list()
	cnt = 1
	for(i in 1:length(partitionDataOut)) {
		if (partitionDataOut[[i]]$n>0) {
			res[[cnt]] = partitionDataOut[[i]]
		cnt = cnt +1
		}	
	}
	return(res)
}

GetConds<-function(tree, catIdxVect, catValueList, idxVarsToExclude = Inf){ 
	# extract all splitting rules from a tree
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
	conds<-list()
	#start by  terminal nodes and recursively trace back
	id.leaves<-which(tree$status==-1) #terminal node's status is -1
	j<-0
	for(i in id.leaves){ #do this for each of the terminal nodes
		j<-j+1
		conds[[j]] = list()
		prevConds<-PrevCond(tree,i, catIdxVect, catValueList)
		
		if (!(prevConds$splitVar %in% idxVarsToExclude)) {
		
			conds[[j]]$rule<-prevConds$cond 
			conds[[j]]$splitVar =  prevConds$splitVar
			conds[[j]]$splitVal = prevConds$splitVal
		} else{
			print("bug: this will give duplicates rules! maybe remove duplicates with hashes at the end; consider using PrunePart")
		}
		
		if(prevConds$motherRow==1) {

			conds[[j]]$pred = tree$prediction[i]
		}
		
		while(prevConds$motherRow>1){ #recur backtracing
			prevConds<-PrevCond(tree,prevConds$motherRow, catIdxVect, catValueList)
			#conds[[j]]$rule<-paste(conds[[j]]$rule," & ",prevConds$cond) #these were commented out because it showed rules from the bottom to the top of the tree
			#conds[[j]]$splitVar = c(conds[[j]]$splitVar , prevConds$splitVar)
			#conds[[j]]$splitVal = c(conds[[j]]$splitVal , prevConds$splitVal)
			if (!(prevConds$splitVar %in% idxVarsToExclude)) {
			
				conds[[j]]$rule<-paste(prevConds$cond, ifelse(is.null(conds[[j]]$rule), "","*") , conds[[j]]$rule) #these replace the three lines above
				conds[[j]]$splitVar = c(prevConds$splitVar, conds[[j]]$splitVar)
				conds[[j]]$splitVal = c(prevConds$splitVal, conds[[j]]$splitVal)
			}
			
			if(prevConds$motherRow==1){
				#conds[[j]]$rule<-paste(conds[[j]]$rule," => ",tree$prediction[i]) #if last recursion reached root node, then add prediction to the sequence of splitting rules
				conds[[j]]$pred = tree$prediction[i]
				break()
			}
		}

	}
	return(conds)
}


PartitionData <- function(getCondsOutORextractCombinedRulesORPrunePartOut, dat) {
	#Partitions dat into disjoint subsets as per the splitting rules in getCondsOut, or extractCombinedRulesOut, or prunePartOut with nonDirectional==F
	#it asssumes the columns in dat correspond to X1, X2, ..., Xp, and that such columns match the order of the predictors used 
	# to train the RF used as input for the RF (that is, match the variable naming in getCondsOutORextractCombinedRulesORPrunePartOut's rules)
	
	input = getCondsOutORextractCombinedRulesORPrunePartOut
	partitionedDat = inputReduced = list()
	cnt = 1
	for (i in 1:length(input)) {
		if (!is.null(input[[i]]$rule)) # we want to remove components of input that are null because we used idxVarsToExclude in getConds
			inputReduced[[cnt]] = input[[i]]
			cnt = cnt +1
	}
	
	for(i in 1:ncol(dat))
		assign(paste("X", i, sep=""), dat[,i]) # note that this will give the output variable column an "X" name, but
												#that is no problem because such "X" will never show up in any rule
	for(i in 1:length(inputReduced)) {
		partitionedDat[[i]] = list()
		partitionedDat[[i]]$datPart = VectToMat(dat[eval(parse(text=inputReduced[[i]]$rule)),])
		partitionedDat[[i]]$rule = inputReduced[[i]]$rule
		#partitionedDat[[i]]$cor =  cor(partitionedDat[[i]]$datPart)
		partitionedDat[[i]]$n = nrow(partitionedDat[[i]]$datPart)
	}	
	
	return(partitionedDat)
}

PrunePart=function(getCondsOut, depth=Inf, nonDirectional = T) {
	# extract *splitting values/variables/rules* up to certain rule depth (number of recursions); if depth ==Inf, it just returns all rules
	# if nonDirectional == T xj>1 and xj<1 are the considered the same rule (i.e., what counts is the splitting value for a given sequence of variable, not the direction)
	# this is to follow Strobl's way of defining partitions
	#ISSUE: a child rule and a parent rule are considered separated (the parent rule is one rule, the child rule is another rule) if nonDirectional = T 
	
	prunedList = prunedListTmp = NULL
	hashVect = NULL
	cnt = 1
	for (i in 1:length(getCondsOut)) {
		prunedListTmp[[i]]=list()
		currLength = length(getCondsOut[[i]]$splitVar)
		prunedListTmp[[i]]$splitVar = getCondsOut[[i]]$splitVar[1:min(currLength,depth)]
		prunedListTmp[[i]]$splitVal = getCondsOut[[i]]$splitVal[1:min(currLength,depth)]


		if(!nonDirectional) # if nonDirectional ==F, we can keep the pruned rules, without having to go through ExtractPart and ExtractCombinedRules
			prunedListTmp[[i]]$rule = paste("as.logical(", paste(strsplit(getCondsOut[[i]]$rule, "*", fixed=T)[[1]][1:min(currLength,depth)], 
				collapse="& "), ")", sep="") 

		if (nonDirectional | depth<Inf){
			#prunedListTmp[[i]]$hash = hashVect[cnt] = paste(c(round(prunedListTmp[[i]]$splitVal,8), prunedListTmp[[i]]$splitVar), sep="", collapse="")
			prunedListTmp[[i]]$hash = hashVect[cnt] = prunedListTmp[[i]]$rule 
			
		}		
		
		currLength = min(currLength,depth)
		prunedListTmp[[i]]$length = currLength
		cnt = cnt + 1

	}
	
	if(nonDirectional | depth<Inf){
		unHashVect = unique(hashVect)
		for (i in 1:length(unHashVect)) { 
			quale = which(hashVect==unHashVect[i])[1]
			prunedList[[i]]=list()
			prunedList[[i]]$splitVar = prunedListTmp[[quale]]$splitVar
			prunedList[[i]]$splitVal = prunedListTmp[[quale]]$splitVal
			prunedList[[i]]$length = prunedListTmp[[quale]]$length
			if(!nonDirectional)
				prunedList[[i]]$rule = prunedListTmp[[quale]]$rule			
		}

		return(prunedList)
	} else{
		return(prunedListTmp)
	}
	
}

ExtractPart = function(prunePartOut, totVarNum) {
	# Returns the union of *splitting values* coming from all rules involving all variables 
	#totVarNum should be = p
	#NOTE: this works for fully numerical instances only!

	vals = list()
	for (i in 1:totVarNum)
		vals[[i]]  = Inf
		
		
	for (i in 1:length(prunePartOut)) {
		for (j in 1:length(prunePartOut[[i]]$splitVar)) {
			whichVar = prunePartOut[[i]]$splitVar[j]
			whichVal = prunePartOut[[i]]$splitVal[j]
			vals[[whichVar]] = c(vals[[whichVar]], whichVal) 
		}
	}
	
	for (i in 1:totVarNum) {
		vals[[i]] = vals[[i]][-1] #removes the Inf, which is used above for coding convenience #### BUGBUG this is ok only for numerical variables
		vals[[i]] = sort(unique(vals[[i]]))
		#if (length(vals[[i]])==0)
		#	vals[[i]] = NULL
	}	
	return(vals)
}

ExtractCombinedRules = function(extractPartOut, condVarVect){
	#Returns all rules obtained via cartesian product of all splitting values of all variables in condVarVect
	#condVarVect: vector of indices corresponding to variables to create partitioning on
	#NOTE: this works for fully numerical instances only!
	rulesList=reducedRulesList=list()
	cnt = 1
	for (i in 1:length(extractPartOut)) {
		if(length(extractPartOut[[i]])>0 &i %in% condVarVect){
			currVar =  paste("X", i, sep="")
			currVals = extractPartOut[[i]]
			currVals = c(-Inf, currVals, Inf)
			rules = NULL
			for (j in 2:(length(currVals))) {
				rules= c(rules, paste(currVar , ">=(", currVals[j-1], ")&", currVar, "<(", currVals[j], ")", sep=""))
			}
			rulesList[[i]] = rules
			if(!is.null(rules)) {
				reducedRulesList[[cnt]] = rules
				cnt = cnt+1
			}
		
		}
	}
	
	if (length(reducedRulesList)>0) { # if none of the variables in condVarVect 
		allCombRules = expand.grid(reducedRulesList)
		allCombRules = apply(allCombRules, 1, paste, collapse="&")
		res= list()

		for(i in 1:length(allCombRules)){
			res[[i]]=list()
			res[[i]]$rule = allCombRules[i]
		}
	} else{
		res = NULL
	}
	return(res)
}




### Functions to transform a random forest object into something human-readable

VectToMat=function(dat, byRow=T) {
	output = dat 
	if(is.vector(dat) ){
		output=as.matrix(dat)
		output=t(output)
	}
	return(output)
}

PropSubTree = function(tree, maxVarIdx, topFrac) {
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
	# returns relative frequencies with which variables are selected as splitting variables ROUGHLY in the top topFrac fraction of splits
	#maxVarIdx should be p (input space dimension)
	splitVarVect = tree[,3] #list of splitting variables (note that the position of a variable in this vector is only a loose indicator of 
							#the level of depth of the tree in which such variable appears a splitting variable)
	splitVarVect = splitVarVect[splitVarVect!=0]
	numSplits = length(splitVarVect)


	splitVarSubVect = splitVarVect[1:round(topFrac*numSplits)]
	tab = table(splitVarSubVect)
	colTab = names(tab)
	tabNum = rep(0,maxVarIdx)
	for (i in 1:length(tab))
		tabNum[as.numeric(colTab[i])] = as.numeric(tab[i])

	tabNum = tabNum/sum(tabNum)

	return(tabNum)
}

PropSubTreeMean = function(rf, topFrac, maxVarIdx) {
	# same as PropSubTree, just averaged out across all trees of a random forest
	#maxVarIdx should be p (input space dimension)
	nTree = rf$ntree
	propMat = NULL
	for(i in 1:nTree){
		tree = getTree(rf,i)
		propMat = rbind(propMat, PropSubTree(tree, maxVarIdx, topFrac))
	}
	res = apply(propMat, 2, mean)
	return(res)
}

PrevCond = function(tree,i, catIdxVect, catValueList){
	
	# trace back splitting rules by one step from node labeled i of tree tree (used by getConds)
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
	# catIdxVect is a vector of indices of categorical input variables in the data
	# catValueList is the a list with p elements  
		# with the j-th element containing a vector of possible levels (relevant only for j %in% catIdxVect) 
	 
	if(i %in% tree[,'right daughter']){
		motherRow<-which(tree[,'right daughter']==i) #the node whose child is node i
		whichVar = tree[motherRow, 'split var']
		whichVal = tree[motherRow, 'split point']
		if(whichVar %in% catIdxVect) {
			whichValIdx = which(as.integer(intToBits(whichVal))==1)
			cond = NULL
			for(k in 1:length(whichValIdx)){ 
				cond= c(cond, paste(paste("X", whichVar, sep=""),"!=", paste("'", catValueList[[whichVar]][whichValIdx[k]], "'", sep=""))) 
			}
			cond = paste(cond, collapse=" & ") ### note that here we want neither variable, as opposed to below, where we'd want any
			
		} else{
			cond<-paste(paste("X", whichVar, sep=""),">=",round(whichVal,7)) #rule that was added in the last recursion to define  node i
		}
		cond = paste("(", cond, ")", sep="")

	}	### if categorical --> category with weight 1 go here, others go right
	
	if(i %in% tree[,'left daughter']){
		motherRow<-which(tree[,'left daughter']==i)
		whichVar = tree[motherRow, 'split var']
		whichVal = tree[motherRow, 'split point']
		if(whichVar %in% catIdxVect) {
			whichValIdx = which(as.integer(intToBits(whichVal))==1)
			cond = NULL
			for(k in 1:length(whichValIdx)){ 
				cond= c(cond, paste(paste("X", whichVar, sep=""),"==", paste("'", catValueList[[whichVar]][whichValIdx[k]], "'", sep=""))) 
			}
			cond = paste(cond, collapse=" | ") 
		} else{
			cond<-paste(paste("X", whichVar, sep=""),"<",round(whichVal,7)) #rule that was added in the last recursion to define  node i
		}
		cond = paste("(", cond, ")", sep="")
	}

	return(list(cond=cond,motherRow=motherRow, splitVar = whichVar, splitVal =whichVal))
}

collapse<-function(x){
	#remove spaces in a word
	
  x<-sub(" ","_",x)

  return(x)
}

ReverseSplit <- function(getCondsOut) { 
	#useless now (was used when getConds was returning splitting rules in inverted order)
	for (i in 1:length(getCondsOut)) {
		getCondsOut[[i]]$splitVar = rev(getCondsOut[[i]]$splitVar)
		getCondsOut[[i]]$splitVal = rev(getCondsOut[[i]]$splitVal)	
	}
	return(getCondsOut)
}

### Older version of functions

PrevCondOld <- function(tree,i){
	# trace back splitting rules by one step from node labeled i of tree tree (used by getConds)
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
  if(i %in% tree[,'right daughter']){
		motherRow<-which(tree[,'right daughter']==i) #the node whose child is node i
		cond<-paste(paste("X", tree[motherRow, 'split var'] , sep=""),">=",round(tree[motherRow, 'split point'],7)) #rule that was added in the last recursion to define  node i
	  }	### if categorical --> category with weight 1 go here, others go right
	  if(i %in% tree[,'left daughter']){
    motherRow<-which(tree[,'left daughter']==i)
		cond<-paste(paste("X", tree[motherRow, 'split var'], sep=""),"<",round(tree[motherRow, 'split point'],7)) #rule that was added in the last recursion to define  node i
  }

  return(list(cond=cond,motherRow=motherRow, splitVar = as.numeric(tree[motherRow, 'split var']), splitVal =tree[motherRow, 'split point']))
}

GetCondsOld<-function(tree){ 
	# extract all splitting rules from a tree
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
	#idxVarsToExlude is a vector of indices of variables for which we want to disregard the (sub) splitting rules

	
	conds<-list()
	#start by  terminal nodes and recursively trace back
	id.leaves<-which(tree$status==-1) #terminal node's status is -1
	j<-0
	for(i in id.leaves){ #do this for each of the terminal nodes
		j<-j+1
		conds[[j]] = list()
		prevConds<-PrevCond(tree,i)
		conds[[j]]$rule<-prevConds$cond 
		conds[[j]]$splitVar =  prevConds$splitVar
		conds[[j]]$splitVal = prevConds$splitVal
		if(prevConds$motherRow==1) {
			# conds[[j]]$rule<-paste(conds[[j]]$rule," => ",tree$prediction[i]) #if this one recursion reached root node, then add prediction to the splitting rule
			conds[[j]]$pred = tree$prediction[i]
		}

		while(prevConds$motherRow>1){ #recur backtracing
			prevConds<-PrevCond(tree,prevConds$motherRow)
			#conds[[j]]$rule<-paste(conds[[j]]$rule," & ",prevConds$cond) #these were commented out because it showed rules from the bottom to the top of the tree
			#conds[[j]]$splitVar = c(conds[[j]]$splitVar , prevConds$splitVar)
			#conds[[j]]$splitVal = c(conds[[j]]$splitVal , prevConds$splitVal)
			conds[[j]]$rule<-paste(prevConds$cond, "&" , conds[[j]]$rule) #these replace the three lines above
			conds[[j]]$splitVar = c(prevConds$splitVar, conds[[j]]$splitVar)
			conds[[j]]$splitVal = c(prevConds$splitVal, conds[[j]]$splitVal)

			if(prevConds$motherRow==1){
				#conds[[j]]$rule<-paste(conds[[j]]$rule," => ",tree$prediction[i]) #if last recursion reached root node, then add prediction to the sequence of splitting rules
				conds[[j]]$pred = tree$prediction[i]
				break()
			}
		}

	}
	return(conds)
}

GetCondsOld2<-function(tree, idxVarsToExclude = Inf){ 
	# extract all splitting rules from a tree
	# tree e.g., as in as.data.frame(getTree(RandomForestObject,1))
	conds<-list()
	#start by  terminal nodes and recursively trace back
	id.leaves<-which(tree$status==-1) #terminal node's status is -1
	j<-0
	for(i in id.leaves){ #do this for each of the terminal nodes
		j<-j+1
		conds[[j]] = list()
		prevConds<-PrevCond(tree,i)
		
		if (!(prevConds$splitVar %in% idxVarsToExclude)) {
		
			conds[[j]]$rule<-prevConds$cond 
			conds[[j]]$splitVar =  prevConds$splitVar
			conds[[j]]$splitVal = prevConds$splitVal
		} else{
			print("bug: this will give duplicates rules! maybe remove duplicates with hashes at the end; consider using PrunePart")
		}
		
		if(prevConds$motherRow==1) {
			# conds[[j]]$rule<-paste(conds[[j]]$rule," => ",tree$prediction[i]) #if this one recursion reached root node, then add prediction to the splitting rule
			conds[[j]]$pred = tree$prediction[i]
		}
		
		while(prevConds$motherRow>1){ #recur backtracing
			prevConds<-PrevCond(tree,prevConds$motherRow)
			#conds[[j]]$rule<-paste(conds[[j]]$rule," & ",prevConds$cond) #these were commented out because it showed rules from the bottom to the top of the tree
			#conds[[j]]$splitVar = c(conds[[j]]$splitVar , prevConds$splitVar)
			#conds[[j]]$splitVal = c(conds[[j]]$splitVal , prevConds$splitVal)
			if (!(prevConds$splitVar %in% idxVarsToExclude)) {
			
				conds[[j]]$rule<-paste(prevConds$cond, ifelse(is.null(conds[[j]]$rule), "","&") , conds[[j]]$rule) #these replace the three lines above
				conds[[j]]$splitVar = c(prevConds$splitVar, conds[[j]]$splitVar)
				conds[[j]]$splitVal = c(prevConds$splitVal, conds[[j]]$splitVal)
			}
			
			if(prevConds$motherRow==1){
				#conds[[j]]$rule<-paste(conds[[j]]$rule," => ",tree$prediction[i]) #if last recursion reached root node, then add prediction to the sequence of splitting rules
				conds[[j]]$pred = tree$prediction[i]
				break()
			}
		}

	}
	return(conds)
}

PrunePartOld <- function(getCondsOut, depth=Inf, nonDirectional = T) {
	# extract *splitting values/variables/rules* up to certain rule depth (number of recursions); if depth ==Inf, it just returns all rules
	# if nonDirectional == T xj>1 and xj<1 are the considered the same rule (i.e., what counts is the splitting value for a given sequence of variable, not the direction)
	# this is to follow Strobl's way of defining partitions
	#ISSUE: a child rule and a parent rule are considered separated (the parent rule is one rule, the child rule is another rule) if nonDirectional = T 
	
	prunedList = prunedListTmp = NULL
	hashVect = NULL
	cnt = 1
	for (i in 1:length(getCondsOut)) {
		prunedListTmp[[i]]=list()
		currLength = length(getCondsOut[[i]]$splitVar)
		prunedListTmp[[i]]$splitVar = getCondsOut[[i]]$splitVar[1:min(currLength,depth)]
		prunedListTmp[[i]]$splitVal = getCondsOut[[i]]$splitVal[1:min(currLength,depth)]

		if (nonDirectional | depth<Inf)
			prunedListTmp[[i]]$hash = hashVect[cnt] = paste(c(round(prunedListTmp[[i]]$splitVal,8), prunedListTmp[[i]]$splitVar), sep="", collapse="")
		
		if(!nonDirectional) # if nonDirectional ==F, we can keep the pruned rules, without having to go through ExtractPart and ExtractCombinedRules
			prunedListTmp[[i]]$rule = paste(strsplit(getCondsOut[[i]]$rule, "&")[[1]][1:min(currLength,depth)], collapse="& ")
			
		prunedListTmp[[i]]$length = currLength
		cnt = cnt + 1

	}
	
	if(nonDirectional | depth<Inf){
		unHashVect = unique(hashVect)
		for (i in 1:length(unHashVect)) { 
			quale = which(hashVect==unHashVect[i])[1] # problem here -- see PrunePart new version
			prunedList[[i]]=list()
			prunedList[[i]]$splitVar = prunedListTmp[[quale]]$splitVar
			prunedList[[i]]$splitVal = prunedListTmp[[quale]]$splitVal
			prunedList[[i]]$length = prunedListTmp[[quale]]$length			
		}
		return(prunedList)
	} else{
		return(prunedListTmp)
	}
	
}

