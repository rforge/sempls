### Armin (first: 2012-09-21; Last: 2013-01-23)

### ic - information criteria
ic <- function(object, LV, criteria, ...){
  UseMethod("ic", object)
}

ic.sempls <- function(object, LV, criteria, ...){
  criteriaOpt <- c("FPE", "AdjRsq", "AIC", "AICu", "AICc", "BIC",
                   "HQ", "HQc", "AIC2", "BIC2", "Cp", "GM")
  if(missing(criteria)) criteria <- criteriaOpt
  stopifnot(criteria %in% criteriaOpt)
  N <- object$N
  pred <- semPLS:::predecessors(object$model)[[LV]]
  pk <- length(pred)
  SSerrk <- semPLS:::deviance.sempls(object, LV)
  sapply(criteria, function(criteria)
         do.call(criteria, list(N, pk, SSerrk, object, LV)))
}

ic_exhaustive <- function(object, LV, criteria, ...){
  UseMethod("ic_exhaustive", object)
}

ic_exhaustive.sempls <- function(object, LV, criteria, ...){
  object <- fitexhaustive(object, LV)
  criteriaOpt <- c("FPE", "AdjRsq", "AIC", "AICu", "AICc", "BIC",
                   "HQ", "HQc", "AIC2", "BIC2", "Cp", "GM")
  if(missing(criteria)) criteria <- criteriaOpt
  stopifnot(criteria %in% criteriaOpt)
  crit <- t(sapply(object, ic, LV = LV, criteria = criteria))
  ret <- list(ic = crit, all_fits = object)
  class(ret) <- "ic_exhaustive"
  return(ret)
}

print.ic_exhaustive <- function(x, minlength = 1, digits = 5, ...){
  rn <- abbreviate(rownames(x$ic), minlength = minlength)
  rownames(x$ic) <- rn
  l <- length(rn)
  indx <- seq_len(log(l + 1, base = 2))
  print(x$ic, digits = digits, ...)
  cat("\n", paste(rn[indx], ": ", names(rn[indx]), "\n", sep = ""), "\n", sep = "")
}

## Final Prediction Error
FPE <- function (N, pk, SSerrk, object, LV){
  SSerrk/N*(1 + (2*pk)/(N-pk))
}


## Adjusted R-sq; This value should match the one given by sempls package.
AdjRsq <- function (N, pk, SSerrk, object, LV){
  SStotal <- N-1
  1 - ((N-1)/(N-pk))*(SSerrk/SStotal)
}

## AIC
AIC <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk/N))
}

AIC2 <- function (N, pk, SSerrk, object, LV){
  ll <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(SSerrk)))
  -2 * ll + (2 * (pk + 1))
}


## Unbiased AIC
AICu <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/(N-pk)) + (2*pk/N))
}


## Corrected AIC
AICc <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (N+pk)/(N-pk-2))
}


## Bayesian IC
BIC <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + pk*log(N)/N)
}

BIC2 <- function (N, pk, SSerrk, object, LV){
  ll <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(SSerrk)))
  -2 * ll + (pk + 1) * log(N)
}


## Hannan Quinn
HQ <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk*log(log(N)))/N)
}

## Corrected Hannan Quinn
HQc <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk*log(log(N)))/(N-pk-2))
}


### need 'saturated' (all main effects) model:

## Mallows Cp
Cp <- function (N, pk, SSerrk, object, LV){
  objectFull <- fitallmain(object, LV)
  MSerrFullModel <-  (1-rSquared(objectFull)[LV,]) # R-Squared from the 'saturated' model.
  (SSerrk / MSerrFullModel) - (N - 2 * pk)
}

## Geweke Meese Criterion
GM <- function (N, pk, SSerrk, object, LV){
  objectFull <- fitallmain(object, LV)
  MSerrFullModel <-  (1-rSquared(objectFull)[LV,]) # R-Squared from the 'saturated' model. 
  (SSerrk/MSerrFullModel) + pk*log(N);
}
