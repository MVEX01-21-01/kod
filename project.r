library(spatstat)

data_moderate <- readRDS('ENDPOINT_DATA/CALF_MODERATE')
data_normal   <- readRDS('ENDPOINT_DATA/CALF_NORMAL')
data <- hyperframe(
  g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
  ppp = c(data_moderate, data_normal)
)

Kbar <- function(ppp, groups, ...) {
  K <- anylapply(ppp, function(ppp) Kest(ppp, ..., ratio=T))
  anylapply(split(K, groups), pool)
}

# Estimating bar(K)
K <- Kbar(data$ppp, data$g, correction='iso')

# Plot and test if they are different
plot(K, sqrt(cbind(pooltheo,pooliso,hiiso,loiso)/pi)-r~r,
     shade=c('hiiso', 'loiso'), equal.scales=T)
print(studpermu.test(data, ppp ~ g))

# Function for fitting cluster processes to K estimate from replicates
repcluster.estK <- function(X, groups, cluster, startpars=NULL, ...) {
  info <- spatstatClusterModelInfo(cluster)
  estK <- match.fun(paste(tolower(cluster), '.estK', sep=''))
  if (!is.null(startpars))
    startpars <- anylapply(startpar, info$checkpar)
  else {
    # just take an average of suggested starting params
    startpars <- anylapply(split(X, groups),
                           function(lppp) rowMeans(sapply(lppp, info$selfstart)))
  }
  
  K <- Kbar(X, groups, correction='best')
  as.anylist(apply(cbind(Kbar=K, par=startpars), 1,
                   function(row) estK(row$Kbar, startpar=row$par, ...)))
}

# Actually fit both models
fit.thomas   <- repcluster.estK(data$ppp, data$g, 'Thomas')
fit.matclust <- repcluster.estK(data$ppp, data$g, 'MatClust')

# spatstat has difficulties computing the bound, let's help
ylim.minconfit.K <- function(fits, cols) {
  m <- sapply(fits, function(fit) {
    ys <- fit$fit[,cols,drop=T]
    c(min(ys, na.rm=T), max(ys, na.rm=T))
  })
  c(min(m[1,]), max(m[2,]))
}
ylim.minconfit.L <- function(fits, cols) {
  m <- sapply(fits, function(fit) {
    ys <- sqrt(fit$fit[,cols,drop=T]/pi) - fit$fit$r
    c(min(ys, na.rm=T), max(ys, na.rm=T))
  })
  c(min(m[1,]), max(m[2,]))
}

# convenience for plotting both
plotfit <- function(f) plot(f, sqrt(cbind(pooltheo,pooliso,fit)/pi)-r~r, ylim=ylim.minconfit.L(f, c('pooliso','pooltheo','fit')))

plotfit(fit.thomas)
plotfit(fit.matclust)
