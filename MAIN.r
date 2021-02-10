# Note: check working directory
library(ggplot2)
theme_set(theme_minimal())
library(spatstat)
source('repcluster.est.r')
source('repcluster.params.r')


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
K <- Kbar(data$ppp, data$g, correction='iso')

# Plot
plot(K, sqrt(cbind(pooltheo,pooliso,hiiso,loiso)/pi)-r~r,
     shade=c('hiiso', 'loiso'), equal.scales=T)
# Test if they are different (not interesting)
#print(studpermu.test(data, ppp ~ g))


# Fitting to aggregate estimate ================================================
fit.thomas   <- repcluster.estK(data$ppp, data$g, 'Thomas')
fit.matclust <- repcluster.estK(data$ppp, data$g, 'MatClust')

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
