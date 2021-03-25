# Estimating bar(lambda) from aggregation recipe in Illian (eq 4.7.1)
intensitybar <- function(X, coeff=1) {
  weighted.mean(sapply(X, intensity)*coeff, sapply(X, area))
}

# Estimating bar(K) from list of patterns
# Standard (ratio of sums weighted by n squared)
Kbar <- function(X, ...) pool(anylapply(X, Kest, ..., ratio=T))
# Linear (ratio of sums weighted by n)
Kbar1 <- function(X, ...) pool(anylapply(X, Kest, ..., ratio=T), weights=sapply(X, npoints))

# Estimating bar(L) from list of patterns
Lbar <- function(X, ...) pool(anylapply(X, Lest, ..., ratio=T))

# Estimating bar(rho) from list of patterns
pcfbar <- function(X, ...) pool(anylapply(X, pcf, ..., ratio=T))

# Fit cluster model to data hyperframe based on estimate
repcluster.estK <- function(hX, cluster, startpars=NULL, lambda=NULL, ...) {
  repcluster_(hX, cluster, 'estK', Kbar, startpars, lambda, ...)
}
repcluster.estpcf <- function(hX, cluster, startpars=NULL, lambda=NULL, ...) {
  repcluster_(hX, cluster, 'estpcf', pcfbar, startpars, lambda, ...)
}
repcluster_ <- function(hX, cluster, estfn, fitfn, startpars=NULL, lambda=NULL, ...) {
  info <- spatstatClusterModelInfo(cluster)
  estfn <- match.fun(paste(tolower(cluster), '.', estfn, sep=''))
  if (!is.null(startpars))
    startpars <- anylapply(startpars, info$checkpar)
  else {
    # just take an average of suggested starting params
    # (as long as we're not too far off, it's ok. This agrees with what
    #  kppm computes for single pattern estimates)
    startpars <- grouped(function(lppp) rowMeans(sapply(lppp, info$selfstart)), hX)
  }
  if (is.null(lambda))
    lambda <- grouped(intensitybar, hX)
  
  S <- grouped(fitfn, hX, correction='best')
  as.anylist(apply(cbind(S=S, par=startpars, lambda=lambda), 1,
                   function(row) estfn(row$S, startpar=row$par, lambda=row$lambda, ...)))
}

# Compute an alternative starting parameter set from branching obs
startpars.branching <- function(dX, pX) {
  c(kappa=intensitybar(pX), scale=2*mean(unlist(sapply(dX, nndist))))
}
