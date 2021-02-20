# Estimating bar(lambda) from aggregation recipe in Illian (eq 4.7.1)
intensitybar <- function(X) {
  intensities <- sapply(X, intensity)
  areas <- sapply(X, area)
  areas.sum <- sum(areas)
  sum(intensities * areas / areas.sum) 
}

# Estimating bar(K) from list of patterns
Kbar <- function(X, ...) pool(anylapply(X, Kest, ..., ratio=T))

# Fit cluster model to data hyperframe based on bar(K) estimate
repcluster.estK <- function(hX, cluster, startpars=NULL, lambda=NULL, ...) {
  info <- spatstatClusterModelInfo(cluster)
  estK <- match.fun(paste(tolower(cluster), '.estK', sep=''))
  if (!is.null(startpars))
    startpars <- anylapply(startpar, info$checkpar)
  else {
    # just take an average of suggested starting params
    # (as long as we're not too far off, it's ok. This agrees with what
    #  kppm computes for single pattern estimates)
    startpars <- grouped(function(lppp) rowMeans(sapply(lppp, info$selfstart)), hX)
  }
  if (is.null(lambda))
    lambda <- grouped(intensitybar, hX)
  
  K <- grouped(Kbar, hX, correction='best')
  as.anylist(apply(cbind(Kbar=K, par=startpars, lambda=lambda), 1,
                   function(row) estK(row$Kbar, startpar=row$par, lambda=row$lambda, ...)))
}