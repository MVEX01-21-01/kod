# Estimating bar(lambda) from aggregation recipe in Illian (eq 4.7.1)
intensitybar <- function(X, groups) {
  anylapply(split(X,groups), function(lppp) {
    intensities <- sapply(lppp, intensity)
    areas <- sapply(lppp, area)
    areas.sum <- sum(areas)
    sum(intensities * areas / areas.sum) 
  })
}

# Estimating bar(K) from data hyperframe based on groups
Kbar <- function(ppp, groups, ...) {
  K <- anylapply(ppp, function(ppp) Kest(ppp, ..., ratio=T))
  anylapply(split(K, groups), pool)
}

# Fit cluster model to data hyperframe based on bar(K) estimate
repcluster.estK <- function(X, groups, cluster, startpars=NULL, lambda=NULL, ...) {
  info <- spatstatClusterModelInfo(cluster)
  estK <- match.fun(paste(tolower(cluster), '.estK', sep=''))
  if (!is.null(startpars))
    startpars <- anylapply(startpar, info$checkpar)
  else {
    # just take an average of suggested starting params
    startpars <- anylapply(split(X, groups),
                           function(lppp) rowMeans(sapply(lppp, info$selfstart)))
  }
  if (is.null(lambda))
    lambda <- intensitybar(X, groups)
  
  K <- Kbar(X, groups, correction='best')
  as.anylist(apply(cbind(Kbar=K, par=startpars, lambda=lambda), 1,
                   function(row) estK(row$Kbar, startpar=row$par, lambda=row$lambda, ...)))
}
