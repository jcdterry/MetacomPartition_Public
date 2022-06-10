## simplified R2 function, that only does adjusted nakagawa and Tjur's D

## always calcualtes speratealty for each species across all sites for the existing data 

## (i.e newdata = NULL, adjust = TRUE , indSite = FALSE,  averageSp = FALSE)

## Also assumed models are probit : will test and stop if not


R2_simplified <- function( hmsc, type = c("nakagawa","tjur" )[1], adjust=TRUE){
  
  Y <- hmsc$data$Y
  
  Xvariance <- apply(hmsc$data$X, 2, var)
  
  if (!any(Xvariance == 0)) { warning("A model without intercept may lead to a wrong R2")}
  if (class(hmsc)[2] != "probit"){stop('This streamlined R2 only does probit models. Use original HMSC::variPart() instead')}    
  if (!any(unique(as.vector(Y)) %in% c(0, 1))) {stop("Only expects binary data")}
  
  nsite <- nrow(Y)
  nsp <- ncol(Y)
  
  Ypred_untransformed <- predict(hmsc, type = "link")  # prediction before link transform, for nakagawa 
  Ypred <- predict(hmsc) # prediction on scale of response (i.e. 0-1), for Tjur's D
  
  if (type == "tjur") {
    Y0 <- Y == 0
    Y1 <- Y == 1
    R2 <- numeric()
    for (i in 1:nsp) {
      R2[i] <- mean(Ypred[Y1[, i], i]) - mean(Ypred[Y0[, i], i])
    }
  }
  
  if (type == "nakagawa") {
    varDist <- 1  ## distribution specific variance (see table 2, N+S)
    YMeans <- matrix(colMeans(Ypred_untransformed), nrow = nsite, ncol = nsp, byrow = TRUE)
    varModelSite <- (Ypred_untransformed - YMeans)^2/(nsite - 1)  # variance in predictions
    varModel <- colSums(varModelSite) # this includes the fixed and the random 
    varAdd <- diag(var(Y - Ypred)) # variance of prediction - observed for each species
    varTot <- varModel + varAdd + varDist
    R2 <- varModel/varTot
  }
  
  if(!adjust){return(R2)}
  
  ## Beginning adjusting, following method in Gelman and Pardoe 2006, Section 3.1
  nexp <- 0
  if (any(names(hmsc$data) == "X")) {
    nexp <- nexp + ncol(hmsc$data$X)    # number of 'fixed' effect parameters estimated per species
  }
  if (any(names(hmsc$data) == "Random")) {
    for (i in 1:ncol(hmsc$results$estimation$latent)) {
      nexp <- nexp + max(sapply(hmsc$results$estimation$latent[, i], ncol))  # find maximum number of latent variables fit across draws
    }
  }
  R2_adjust <- 1 - ((nsite - 3)/(nsite - nexp - 2)) * (1 -R2)
  return(R2_adjust)
}
