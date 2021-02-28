# Note: check working directory
library(ggplot2)
theme_set(theme_minimal())
library(spatstat)

# Load data ====================================================================
source('load_data.r')

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

# Analysing aggregate summary statistics =======================================
source('repcluster.r')

# Estimate bar(K) for each group
K <- grouped(Kbar, data, correction='iso')

# Plot
plot(K, sqrt(cbind(pooltheo,pooliso,hiiso,loiso)/pi)-r~r,
     shade=c('hiiso', 'loiso'), equal.scales=T)

# Estimate bar(rho) for each group
pcfb <- grouped(pcfbar, data, correction='iso')

# Plot
plot(pcfb, cbind(pooltheo,pooliso,hiiso,loiso)~r,
     shade=c('hiiso', 'loiso'), equal.scales=T)


# Fitting to aggregate estimate ================================================
fit.thomas.K   <- repcluster.estK(data, 'Thomas')
fit.matclust.K <- repcluster.estK(data, 'MatClust')
fit.thomas.g   <- repcluster.estpcf(data, 'Thomas')
fit.matclust.g <- repcluster.estpcf(data, 'MatClust')

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
plotfit.L <- function(f) plot(f, sqrt(cbind(pooltheo,pooliso,fit)/pi)-r~r,
                            ylim=ylim.minconfit.L(f, c('pooliso','pooltheo','fit')))
plotfit.g <- function(f) plot(f, cbind(pooltheo,pooliso,fit)~r)
plotfit.L(fit.thomas.K)
plotfit.L(fit.matclust.K)
plotfit.g(fit.thomas.g)
plotfit.g(fit.matclust.g)

# Envelope test ----
source('multiGET.r')

# Compute the envelopes
# (Can be done using multiprocessing: check batch_envelopes.r)
handlers(handler_rstudio())
with_progress({
  p <- progressor(2*length(data$ppp))
  rp <- progress_aggregator(p)
  rp(envs.thomas   <- grouped(multiGET.composite, data, fit.thomas.K, c(Gest), alpha=0.1, type='erl'))
  rp(envs.matclust <- grouped(multiGET.composite, data, fit.matclust.K, c(Gest), alpha=0.1, type='erl'))
})
multiGET.plot(envs.thomas)
multiGET.plot(envs.matclust)
