# Note: check working directory
library(ggplot2)
theme_set(theme_minimal())
library(spatstat)

# Helper function splitting over groups.
# This enables us to write functions in terms of a single set of patterns,
# and then using this function compose to groups
#grouped <- function(f, data, ...) anylapply(split(data$ppp, data$g), f, ...)
grouped <- function(f, data, ..., SIMPLIFY=F) {
  as.anylist(mapply(f, split(data$ppp, data$g), ..., SIMPLIFY=SIMPLIFY))
}


# Load data ====================================================================
data_moderate <- readRDS('ENDPOINT_DATA/CALF_MODERATE')
data_normal   <- readRDS('ENDPOINT_DATA/CALF_NORMAL')
data <- hyperframe(
  g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
  ppp = c(data_moderate, data_normal)
)


# Box plot factored by group of per observation parameter fits =================
params.each <- function(X, cluster) {
  sapply(mapply(kppm, X=X, cluster=cluster, SIMPLIFY=F), function(fit) fit$clustpar)
}
fit.each.thomas   <- params.each(data$ppp, 'Thomas')
fit.each.matclust <- params.each(data$ppp, 'MatClust')
longparams <- function(pars, model, group) {
  # reshape the data into easily visualized format
  df <- data.frame(t(pars), model, group)
  varying <- colnames(df)[1:(ncol(df)-2)]
  df <- reshape(df, direction='long', varying=varying, v.names='value',
          timevar='param', times=varying)
  rownames(df) <- NULL
  df
}

data.params <- rbind(longparams(fit.each.thomas, 'Thomas', data$g),
                     longparams(fit.each.matclust, 'MatClust', data$g))
ggplot(data.params, aes(group, value)) +
  geom_boxplot() + facet_wrap(~ model + param, scales='free')


# Estimating bar(K) ============================================================
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

K <- grouped(Kbar, data, correction='iso')

# Plot
plot(K, sqrt(cbind(pooltheo,pooliso,hiiso,loiso)/pi)-r~r,
     shade=c('hiiso', 'loiso'), equal.scales=T)


# Fitting to aggregate estimate ================================================
fit.thomas   <- repcluster.estK(data, 'Thomas')
fit.matclust <- repcluster.estK(data, 'MatClust')

# Plotting the fits ----
# spatstat seems to have difficulties with determining limits
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
plotfit <- function(f) plot(f, sqrt(cbind(pooltheo,pooliso,fit)/pi)-r~r,
                            ylim=ylim.minconfit.L(f, c('pooliso','pooltheo','fit')))
plotfit(fit.thomas)
plotfit(fit.matclust)


# Envelope test ----
source('multiGET.r')

# Compute the envelopes
tic <- proc.time()
envs.thomas   <- grouped(multiGET.composite, data, fit.thomas, c(Gest), alpha=0.05, type='area')
envs.matclust <- grouped(multiGET.composite, data, fit.matclust, c(Gest), alpha=0.05, type='area')
proc.time() - tic
multiGET.plot(envs.thomas)
multiGET.plot(envs.matclust)
