variPart_returnboth <- function (hmsc, groupX, HMSCprior = NULL, type = "I", indSite = FALSE, 
          R2adjust = FALSE, verbose = TRUE, ...) {
  ####
  ## This is identical to HMSC::variPart(), just that it returns both the R2 and the overlap fractions for the type III only (Type 1 deleted for brevity)
  
  if (type != "III") {
    stop("type needs to be  'III'")
  }
  nsite <- nrow(hmsc$data$Y)
  nsp <- ncol(hmsc$data$Y)
  nX <- ncol(hmsc$data$X)
  nAuto <- length(hmsc$data$Auto)
  nRandom <- ncol(hmsc$data$Random)
  if (!is.null(nX)) {
    if (length(groupX) != nX) {
      stop("groupX should have the same length as there are variables in X (including the intercept)")
    }
  }
  if (type == "III") {
    family <- attributes(hmsc)$class[2]
    model <- hmsc
    remove(hmsc)
    groupXNames <- unique(groupX)
    RandomNames <- names(model$data$Random)
    AutoNames <- names(model$data$Auto)
    setsVar <- c(groupXNames, RandomNames, AutoNames)
    nsetsVar <- length(setsVar)
    subModelBase <- vector("list", length = nsetsVar - 1)
    names(subModelBase) <- paste("var", 1:(nsetsVar - 1), 
                                 sep = "")
    for (i in 1:(nsetsVar - 1)) {
      subModelBase[[i]] <- combn(setsVar, i)
    }
    X <- vector("list", length = nsetsVar - 1)
    names(X) <- paste("var", 1:(nsetsVar - 1), sep = "")
    for (i in 1:(nsetsVar - 1)) {
      X[[i]] <- vector("list", length = ncol(subModelBase[[i]]))
    }
    for (i in 1:(nsetsVar - 1)) {
      for (j in 1:ncol(subModelBase[[i]])) {
        if (is.null(model$data$X)) {
          X[[i]][[j]] <- NULL
        }
        else {
          X[[i]][[j]] <- matrix(NA, nrow = nrow(model$data$X), 
                                ncol = 0)
        }
      }
    }
    Random <- vector("list", length = nsetsVar - 1)
    names(Random) <- paste("var", 1:(nsetsVar - 1), sep = "")
    for (i in 1:(nsetsVar - 1)) {
      Random[[i]] <- vector("list", length = ncol(subModelBase[[i]]))
    }
    for (i in 1:(nsetsVar - 1)) {
      for (j in 1:ncol(subModelBase[[i]])) {
        if (is.null(model$data$Random)) {
          Random[[i]][[j]] <- NULL
        }
        else {
          Random[[i]][[j]] <- as.data.frame(matrix(NA, 
                                                   nrow = nrow(model$data$Random), ncol = nRandom))
          colnames(Random[[i]][[j]]) <- RandomNames
        }
      }
    }
    Auto <- vector("list", length = nsetsVar - 1)
    names(Auto) <- paste("var", 1:(nsetsVar - 1), sep = "")
    for (i in 1:(nsetsVar - 1)) {
      Auto[[i]] <- vector("list", length = ncol(subModelBase[[i]]))
    }
    for (i in 1:(nsetsVar - 1)) {
      for (j in 1:ncol(subModelBase[[i]])) {
        Auto[[i]][[j]] <- vector("list", length = nAuto)
      }
    }
    for (i in 1:(nsetsVar - 1)) {
      for (j in 1:ncol(subModelBase[[i]])) {
        if (any(subModelBase[[i]][, j] %in% groupXNames)) {
          nX <- sum(subModelBase[[i]][, j] %in% groupXNames)
          XSel <- which(groupXNames %in% subModelBase[[i]][, 
                                                           j])
          namesX <- character()
          for (k in 1:nX) {
            X[[i]][[j]] <- cbind(X[[i]][[j]], model$data$X[, 
                                                           which(groupX == groupXNames[XSel[k]])])
            namesX <- c(namesX, colnames(model$data$X)[which(groupX == 
                                                               groupXNames[XSel[k]])])
          }
          colnames(X[[i]][[j]]) <- namesX
        }
        else {
          X[[i]][[j]] <- list(NULL)
        }
        if (any(subModelBase[[i]][, j] %in% RandomNames)) {
          nRandom <- sum(subModelBase[[i]][, j] %in% 
                           RandomNames)
          RandomSel <- which(RandomNames %in% subModelBase[[i]][, 
                                                                j])
          for (k in 1:nRandom) {
            Random[[i]][[j]][, k] <- model$data$Random[, 
                                                       which(RandomNames == RandomNames[RandomSel[k]])]
          }
        }
        else {
          Random[[i]][[j]] <- list(NULL)
        }
        if (any(is.na(Random[[i]][[j]]))) {
          colToRm <- unique(which(is.na(Random[[i]][[j]]), 
                                  arr.ind = TRUE)[, 2])
          Random[[i]][[j]] <- Random[[i]][[j]][, -colToRm]
        }
        if (any(subModelBase[[i]][, j] %in% AutoNames)) {
          nAuto <- sum(subModelBase[[i]][, j] %in% AutoNames)
          AutoSel <- which(AutoNames %in% subModelBase[[i]][, 
                                                            j])
          for (k in 1:nAuto) {
            Auto[[i]][[j]][[k]] <- model$data$Auto[[which(AutoNames == 
                                                            AutoNames[AutoSel[k]])]]
            names(Auto[[i]][[j]])[k] <- names(model$data$Auto)[k]
          }
        }
        if (length(Auto[[i]][[j]]) > 1) {
          AutoNull <- which(sapply(Auto[[i]][[j]], is.null))
          if (length(AutoNull) > 0) {
            Auto[[i]][[j]][AutoNull] <- NULL
          }
        }
      }
    }
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
      abs(x - round(x)) < tol
    }
    if (any(names(model$data) == "X")) {
      iterNames <- dimnames(model$results$estimation$paramX)[[3]]
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      burnNames <- dimnames(model$results$burning$paramX)[[3]]
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)
      thinModel <- iterVal[2] - iterVal[1]
    }
    if (any(names(model$data) == "Random")) {
      iterNames <- rownames(model$results$estimation$latent)
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      burnNames <- rownames(model$results$burning$latent)
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)
      thinModel <- iterVal[2] - iterVal[1]
    }
    if (any(names(model$data) == "Auto")) {
      iterNames <- rownames(model$results$estimation$latentAuto)
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      burnNames <- rownames(model$results$burning$latentAuto)
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)
      thinModel <- iterVal[2] - iterVal[1]
    }
    if (!is.wholenumber(niterModel/thinModel)) {
      niterModel <- niterModel + thinModel - 1
    }
    if (!is.wholenumber(nburnModel/thinModel)) {
      nburnModel <- nburnModel + thinModel - 1
    }
    if (!indSite) {
      R2model <- vector("list", length = nsetsVar)
      names(R2model) <- paste("var", 1:nsetsVar, sep = "")
      for (i in 1:(nsetsVar - 1)) {
        R2model[[i]] <- matrix(NA, nrow = nsp, ncol = ncol(subModelBase[[i]]))
        rownames(R2model[[i]]) <- colnames(model$data$Y)
        columnNames <- character()
        for (j in 1:ncol(subModelBase[[i]])) {
          columnNames[j] <- paste(subModelBase[[i]][, 
                                                    j], collapse = "-")
        }
        colnames(R2model[[i]]) <- columnNames
      }
      R2model[[nsetsVar]] <- matrix(NA, nrow = nsp, ncol = 1)
      colnames(R2model[[nsetsVar]]) <- paste(setsVar, collapse = "-")
      rownames(R2model[[nsetsVar]]) <- colnames(model$data$Y)
    }
    else {
      R2model <- vector("list", length = nsetsVar)
      names(R2model) <- paste("var", 1:nsetsVar, sep = "")
      for (i in 1:(nsetsVar - 1)) {
        R2model[[i]] <- array(dim = c(nsite, nsp, ncol(subModelBase[[i]])))
        rownames(R2model[[i]]) <- rownames(model$data$Y)
        colnames(R2model[[i]]) <- colnames(model$data$Y)
        dim3Names <- character()
        for (j in 1:ncol(subModelBase[[i]])) {
          dim3Names[j] <- paste(subModelBase[[i]][, j], 
                                collapse = "-")
        }
        dimnames(R2model[[i]])[[3]] <- dim3Names
      }
      R2model[[nsetsVar]] <- array(dim = c(nsite, nsp, 
                                           1))
      rownames(R2model[[nsetsVar]]) <- rownames(model$data$Y)
      colnames(R2model[[nsetsVar]]) <- colnames(model$data$Y)
      dimnames(R2model[[nsetsVar]])[[3]] <- paste(setsVar, 
                                                  collapse = "-")
    }
    if (verbose) {
      print(paste(sum(sapply(subModelBase, ncol)), "submodels need to be estimated"))
    }
    counter <- 1
    for (i in 1:(nsetsVar - 1)) {
      for (j in 1:ncol(subModelBase[[i]])) {
        if (all(sapply(X[[i]][[j]], is.null))) {
          XUse <- NULL
        }
        else {
          XUse <- X[[i]][[j]]
        }
        if (all(sapply(Random[[i]][[j]], is.null))) {
          RandomUse <- NULL
        }
        else {
          RandomUse <- Random[[i]][[j]]
        }
        if (all(sapply(Auto[[i]][[j]], is.null))) {
          AutoUse <- NULL
        }
        else {
          AutoUse <- Auto[[i]][[j]]
        }
        options(warn = -1)
        HMSCdataObj <- as.HMSCdata(Y = model$data$Y, 
                                   X = XUse, Random = RandomUse, Auto = AutoUse, 
                                   interceptX = TRUE, scaleX = FALSE)
        options(warn = 0)
        submodel <- hmsc(data = HMSCdataObj, priors = HMSCprior, 
                         family = family, niter = niterModel, nburn = nburnModel, 
                         thin = thinModel, verbose = FALSE, ...)
        if (!indSite) {
          R2model[[i]][, j] <- Rsquared(submodel, type = "nakagawa", 
                                        adjust = R2adjust, averageSp = FALSE)
        }
        else {
          R2model[[i]][, , j] <- Rsquared(submodel, type = "nakagawa", 
                                          adjust = R2adjust, indSite = indSite, averageSp = FALSE)
        }
        if (verbose) {
          print(paste("Number of submodels estimated:", 
                      counter))
        }
        counter <- counter + 1
      }
    }
    if (!indSite) {
      R2model[[nsetsVar]][, 1] <- Rsquared(model, type = "nakagawa", 
                                           adjust = R2adjust, indSite = indSite, averageSp = FALSE)
    }
    else {
      R2model[[nsetsVar]][, , 1] <- Rsquared(model, type = "nakagawa", 
                                             adjust = R2adjust, indSite = indSite, averageSp = FALSE)
    }
    if (!indSite) {
      fraction <- vector("list", length = nsetsVar - 1)
      names(fraction) <- paste("overlap", (nsetsVar - 1):1, 
                               sep = "")
      for (i in 1:(nsetsVar - 1)) {
        fraction[[i]] <- matrix(NA, nrow = nsp, ncol = ncol(subModelBase[[i]]))
        rownames(fraction[[i]]) <- colnames(model$data$Y)
      }
    }
    else {
      fraction <- vector("list", length = nsetsVar - 1)
      names(fraction) <- paste("overlap", (nsetsVar - 1):1, 
                               sep = "")
      for (i in 1:(nsetsVar - 1)) {
        fraction[[i]] <- array(dim = c(nsite, nsp, ncol(subModelBase[[i]])))
        rownames(fraction[[i]]) <- rownames(model$data$Y)
        colnames(fraction[[i]]) <- colnames(model$data$Y)
      }
    }
    if (!indSite) {
      ref <- unlist(strsplit(colnames(R2model[[nsetsVar]]), 
                             "-"))
      for (i in 1:(nsetsVar - 1)) {
        columnNames <- character()
        for (j in 1:ncol(subModelBase[[i]])) {
          comp <- unlist(strsplit(colnames(R2model[[i]])[j], 
                                  "-"))
          fraction[[i]][, j] <- R2model[[nsetsVar]] - 
            R2model[[i]][, j]
          columnNames[j] <- paste(ref[which(!(ref %in% 
                                                comp))], collapse = "-")
        }
        colnames(fraction[[i]]) <- columnNames
      }
      if (nsetsVar > 2) {
        for (i in 1:(nsetsVar - 2)) {
          for (j in 1:ncol(fraction[[i]])) {
            fracToSub <- unlist(strsplit(colnames(fraction[[i]])[j], 
                                         "-"))
            fracSel <- colnames(fraction[[nsetsVar - 
                                            1]]) %in% fracToSub
            fracSub <- fraction[[nsetsVar - 1]][, fracSel]
            if (is.null(ncol(fracSub))) {
              fracSub <- matrix(fracSub, nrow = ncol(model$data$Y))
            }
            fraction[[i]][, j] <- fraction[[i]][, j] - 
              rowSums(fracSub)
          }
        }
      }
      overlapAll <- drop(R2model[[nsetsVar]]) - rowSums(sapply(fraction, 
                                                               rowSums))
      overlapAll <- matrix(overlapAll, ncol = 1)
      rownames(overlapAll) <- colnames(model$data$Y)
      colnames(overlapAll) <- paste(setsVar, collapse = "-")
      fraction[[paste("overlap", nsetsVar, sep = "")]] <- overlapAll
      fraction <- fraction[order(names(fraction))]
    }
    if (indSite) {
      ref <- unlist(strsplit(dimnames(R2model[[nsetsVar]])[[3]], 
                             "-"))
      for (i in 1:(nsetsVar - 1)) {
        dim3Names <- character()
        for (j in 1:ncol(subModelBase[[i]])) {
          comp <- unlist(strsplit(dimnames(R2model[[i]])[[3]][j], 
                                  "-"))
          fraction[[i]][, , j] <- R2model[[nsetsVar]][, 
                                                      , 1] - R2model[[i]][, , j]
          dim3Names[j] <- paste(ref[which(!(ref %in% 
                                              comp))], collapse = "-")
        }
        dimnames(fraction[[i]])[[3]] <- dim3Names
      }
      if (nsetsVar > 2) {
        for (i in 1:(nsetsVar - 2)) {
          for (j in 1:dim(fraction[[i]])[[3]]) {
            fracToSub <- unlist(strsplit(dimnames(fraction[[i]])[[3]][j], 
                                         "-"))
            fracSel <- dimnames(fraction[[nsetsVar - 
                                            1]])[[3]] %in% fracToSub
            fracSub <- fraction[[nsetsVar - 1]][, , fracSel]
            if (is.null(ncol(fracSub))) {
              fracSub <- array(fracSub, nrow = nsite, 
                               ncol = nsp)
            }
            fraction[[i]][, , j] <- fraction[[i]][, , 
                                                  j] - apply(fracSub, 1:2, sum)
          }
        }
      }
      fracSum <- lapply(fraction, function(x) apply(x, 
                                                    1:2, sum))
      fracSumArray <- array(NA, dim = c(nsite, nsp, length(fracSum)))
      for (i in 1:length(fracSum)) {
        fracSumArray[, , i] <- fracSum[[i]]
      }
      overlapAll <- drop(R2model[[nsetsVar]]) - apply(fracSumArray, 
                                                      1:2, sum)
      overlapAll <- matrix(overlapAll, nrow = nsite, ncol = nsp)
      rownames(overlapAll) <- rownames(model$data$Y)
      colnames(overlapAll) <- colnames(model$data$Y)
      fraction[[paste("overlap", nsetsVar, sep = "")]] <- overlapAll
      fraction <- fraction[order(names(fraction))]
    }
    res <- fraction
  }
  return(list(res=res, R2model=R2model))
}


