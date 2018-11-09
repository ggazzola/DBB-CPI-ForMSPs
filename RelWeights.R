relweights <- function(fit,...){                         
  R <- cor(fit$model)   
  nvar <- ncol(R)          
  rxx <- R[2:nvar, 2:nvar] 
  rxy <- R[2:nvar, 1]      
  svd <- eigen(rxx)        
  evec <- svd$vectors                           
  ev <- svd$values         
  delta <- diag(sqrt(ev))  
  lambda <- evec %*% delta %*% t(evec)        
  lambdasq <- lambda ^ 2   
  beta <- solve(lambda) %*% rxy           
  rsquare <- colSums(beta ^ 2)                   
  rawwgt <- lambdasq %*% beta ^ 2    
  #import <- (rawwgt / rsquare) * 100 
  import <- (rawwgt / rsquare)  #GGG

  #lbls <- names(fit$model[2:nvar])   
  lbls <- paste("X", 1:(nvar-1), sep=".")   #GGG
	rownames(import) <- lbls
  colnames(import) <- "Weights"
  #barplot(t(import),names.arg=lbls,
   #  ylab="% of R-Square",
    # xlab="Predictor Variables",
     #main="Relative Importance of Predictor Variables", 
     #sub=paste("R-Square=", round(rsquare, digits=3)),
     #...)  
return(import)
}
# using relweights()
#fit <- lm(Murder ~ Population + Illiteracy + Income + 
 #   Frost, data = states)
#relweights(fit, col = "lightgrey")