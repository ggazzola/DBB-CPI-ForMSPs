require(relaimpo)	#package 'relaimpo' should be installed prior to runnning this code
dat=read.csv("Data2.csv")
location = dat[,1]
dat = dat[location<4,-1]
formList = impList = list()
formList[[1]]=as.formula("y1~x1+x2+x28")
formList[[2]]=as.formula("x28~x3+x4+x5+x17+x20+x22+x26")
formList[[3]]=as.formula("x17~x3+x4+x5+x6+x7+x12+x14+x15") #change in graph
formList[[4]]=as.formula("x20~x5")


for (i in 1:4){
	impList[[i]] = list(formula = NULL, model = NULL, rsq = NULL, relImp = NULL)
	impList[[i]]$model = lm(formList[[i]], data = dat)	# estimation of linear model for the i-th stage
	if(i<4){
		 impList[[i]]$relImp = calc.relimp(impList[[i]]$model, type = "lmg", rela=T)@lmg #computation of regressors' relative importance using general dominance weights ("lmg" method in package relaimpo)
	} else{
		impList[[i]]$relImp = 1
	}
	impList[[i]]$formula = formList[[i]]
	impList[[i]]$rsq = summary(impList[[i]]$model)$r.squared
	impList[[i]]$relImp = impList[[i]]$relImp* impList[[i]]$rsq
	impList[[i]]$eps = 1-impList[[i]]$rsq
	
}

importList = list()

importList[[1]] = impList[[1]]$relImp[1]
importList[[2]] = impList[[1]]$relImp[2]
importList[[28]] = impList[[1]]$relImp[3]*impList[[2]]$eps
importList[[3]] = impList[[1]]$relImp[3]*impList[[2]]$eps