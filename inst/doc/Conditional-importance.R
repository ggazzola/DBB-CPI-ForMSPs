### R code from vignette source 'Conditional-importance.rnw'

###################################################
### code chunk number 1: c1
###################################################
options(width=60)
options(warn=-1)


###################################################
### code chunk number 2: c2
###################################################
require(extendedForest)
require(MASS)

#Set up covariance matrix

Cov <- matrix(0,12,12)
Cov[1:4,1:4] <- 0.9
diag(Cov)[] <- 1

#Coefficients for linear model

beta <- c(5,5,2,0,-5,-5,-2,0,0,0,0,0)

# Set the maximum number of partitions to compute the importance 
# from conditional permutation distribution of each variable
maxK<-c(0,2,4,6)

# Set the number of records (or sites) and the number of simulations.
nsites<- 100
nsim <- 100

imp <- array(0,dim=c(12,4,nsim))

#Simulation 

set.seed(222)

for (sim in 1:nsim) { 
  X <- mvrnorm(nsites, rep(0,12), Sigma=Cov)
  Y <- X%*%beta + rnorm(nsites,0,0.5)
  df <- cbind(Y=Y,as.data.frame(X))
  for (lev in 1:4) { 
    RF <- randomForest(Y ~ .,df, maxLevel=maxK[lev], importance=TRUE, ntree=500, corr.threshold=0.5,mtry=8)
    imp[,lev,sim] <- RF$importance[,1]
  }
}
dimnames(imp) <- list(rownames(RF$importance), as.character(maxK), NULL)
imp <- as.data.frame.table(imp)


###################################################
### code chunk number 3: c3
###################################################
require(lattice)
names(imp) <- c("var","maxK","sim","importance")
print(bwplot(var ~ importance | ifelse(maxK=="0", "Marginal", paste("Conditional: Level",maxK,sep="=")),imp,  as.table=T))


###################################################
### code chunk number 4: Conditional-importance.rnw:214-215
###################################################
toLatex(sessionInfo())


