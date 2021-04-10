### MVEX01-21-01
### This file should render the results of the report exactly,
### for reproducibility.

# Initialization ==========
if (!dir.exists('DATA_ENFS')) {
  stop('Can\'t find pattern subdir DATA_ENFS. Check your working directory.')
}

# This is important for reproducibility!
# But note that the pre-stored envelopes also have to be reproducibly generated.
set.seed(012101)

library(ggplot2)
theme_set(theme_minimal())
library(spatstat)

source('util.r')
source('repcluster.r')
source('multiGET.r')
source('plotting.r')
data <- loaddata()
data.branching <- loaddata.branching()

# Clean output directory
unlink('report_out', recursive=T)
dir.create('report_out')

# Helper to unnest single envelopes
unlist1 <- function(l) lapply(l, function(i) i[[1]])

# Helper to write csvsimple-safe csv
write.csvX <- function(...) write.csv(..., row.names=F, quote=F)

# General pattern information
write.csvX(data.frame(
  id=1:length(data$ppp),
  group=data$g,
  n=sapply(data$ppp, npoints),
  xmin=sapply(data$ppp, function(X) X$window$xrange[1]),
  xmax=sapply(data$ppp, function(X) X$window$xrange[2]),
  ymin=sapply(data$ppp, function(X) X$window$yrange[1]),
  ymax=sapply(data$ppp, function(X) X$window$yrange[2]),
  absW=sapply(data$ppp, area)
), 'report_out/info.csv')

# CSR test on patterns ==========
envs.csr <- csrenvs(data$ppp, nsim=999)
g <- plot.envs.single(envs.csr, ncol=8)
ggsave('report_out/01_csr.pdf', plot=g, width=5.5, height=3)

# Individual parameter fits ==========
fit.each.thomas   <- params.each(data$ppp, 'Thomas')
fit.each.matclust <- params.each(data$ppp, 'MatClust')

# Output csv of parameters
write.csvX(data.frame(id=1:nrow(data),thomas=t(fit.each.thomas),matclust=t(fit.each.matclust)), 'report_out/ind.par.csv')

# Helper function to rescale to kappa * area
paramrescale <- function(params, area) {
  params2 <- as.data.frame(params)
  params2['kappa',] <- params['kappa',] * area
  rownames(params2)[rownames(params2) == 'kappa'] <- 'kappa*area'
  params2
}

# Helper function to reshape parameter data into a format better for ggplot
longparams <- function(pars, model, group, area=NULL) {
  if (!is.null(area))
    pars <- paramrescale(pars, area)
  df <- data.frame(t(pars), model, group, check.names=F)
  varying <- colnames(df)[1:(ncol(df)-2)]
  df <- reshape(df, direction='long', varying=varying, v.names='value',
                timevar='param', times=varying)
  rownames(df) <- NULL
  df
}

# Calculate areas and parameter df
areas <- sapply(data$ppp, area)
data.params <- rbind(longparams(fit.each.thomas, 'Thomas', data$g, area=areas),
                     longparams(fit.each.matclust, 'MatClust', data$g, area=areas))
g.indparams <- ggplot(data.params, aes(group, value)) +
  geom_boxplot() + facet_wrap(~ model + param, scales='free')
ggsave('report_out/02_ind.box.pdf', plot=g.indparams, width=5.2, height=5.2)

# Individual envelopes ==========
# Load pre-generated envelope, from batch_single_envelopes.r
# TODO fix data format?
envs.ind <- loadenv('envs199_single.rds')
g <- plot.envs.single(unlist1(envs.ind$Thomas), ncol=8)
ggsave('report_out/03_envs.ind.Thomas.pdf', plot=g, width=5.5, height=3)
g <- plot.envs.single(unlist1(envs.ind$MatClust), ncol=8)
ggsave('report_out/04_envs.ind.MatClust.pdf', plot=g, width=5.5, height=3)

# Group parameter fits ==========
## Show the variance in individual K, and the fitted curves
K <- grouped(Kbar, data, correction='iso')
fit.thomas.K <- repcluster.estK(data, 'Thomas')
fit.matclust.K <- repcluster.estK(data, 'MatClust')
write.csvX(data.frame(
  group=names(fit.thomas.K),
  thomas=t(mapply(function(x) x$modelpar, fit.thomas.K)),
  matclust=t(mapply(function(x) x$modelpar, fit.matclust.K))),
  'report_out/bar.par.csv')

brpars <- Map(startpars.branching, split(data$ppp, data$g), split(data.branching$ppp, data.branching$g))
fit.thomas.K.br <- repcluster.estK(data, 'Thomas', startpars=brpars)
fit.matclust.K.br <- repcluster.estK(data, 'MatClust', startpars=brpars)

Llimits <- list(
  xmin = max(sapply(K, function(fv) min(fv$r))),
  xmax = min(sapply(K, function(fv) max(fv$r))),
  ymin = min(sapply(K, function(fv) min(sqrt(fv$loiso/pi)-fv$r, sqrt(fv$lotheo/pi)-fv$r, na.rm=T))),
  ymax = max(sapply(K, function(fv) max(sqrt(fv$hiiso/pi)-fv$r, sqrt(fv$hitheo/pi)-fv$r, na.rm=T)))
)
Lgs <- Map(function(fv, th, mc) {
  df <- as.data.frame(fv)
  ggplot(df, aes(r, sqrt(pooliso/pi)-r)) +
    geom_line(aes(colour='empirical')) +
    geom_ribbon(aes(ymin=sqrt(loiso/pi)-r,ymax=sqrt(hiiso/pi)-r),alpha=0.1, colour=c("#e0e0e0")) +
    geom_line(aes(,sqrt(pooltheo/pi)-r), linetype = 'dashed') +
    geom_line(aes(,sqrt(th$fit$fit/pi)-r, colour='Thomas'), linetype='dashed') +
    geom_line(aes(,sqrt(mc$fit$fit/pi)-r, colour='MatClust'), linetype='dashed') +
    coord_cartesian(xlim=c(Llimits$xmin, Llimits$xmax), ylim=c(Llimits$ymin, Llimits$ymax)) +
    scale_colour_manual(values=c(c("#444444"),c("#F7951F"),c("#0154A6"))) +
    xlab(expression(italic(r))) + ylab(expression(italic(L(r)-r))) +
    theme(legend.position='bottom', legend.title=element_blank())
}, K, fit.thomas.K, fit.matclust.K)
Laux <- auxgrobs(Lgs[[1]])
g <- combi.grouped(lapply(Lgs, function(g) list(g + labs(x=NULL, y=NULL) + theme(legend.position='none'))), aux=Laux)
ggsave('report_out/05_bar.L.fit.pdf', plot=g, width=5, height=3)

## Overlay the group results on the boxplots
df.fit <- rbind(
  longparams(rbind(sapply(fit.thomas.K, function(f) f$modelpar)), 'Thomas', c('MODERATE','NORMAL'), area=mean(areas)),
  longparams(rbind(sapply(fit.matclust.K, function(f) f$modelpar)), 'MatClust', c('MODERATE','NORMAL'), area=mean(areas))
)
df.fit[df.fit=='sigma' | df.fit=='R'] <- 'scale'
df.true <- data.frame(
  model=rep(c('Thomas','MatClust'),each=4),
  group=rep(c('MODERATE','NORMAL'),each=2,length.out=8),
  param=rep(c('kappa*area','mu'),length.out=8),
  value=rep(unlist(grouped(function(X, dX) {
    ib <- intensitybar(X)
    # NOTE: reweighted to kappa*area
    c(intensitybar(X,coeff=sapply(X, area)), intensitybar(dX,coeff=1/ib))
  }, data.branching, dX=split(data$ppp, data$g))),length.out=8)
)
g <- g.indparams +
  geom_point(data=df.fit, shape=3, color='red') +
  geom_point(data=df.true, shape=4, color='blue')
ggsave('report_out/06_bar.box.pdf', plot=g, width=5.2, height=5.2)

# Group envelopes ==========
# Load pre-generated envelope, from batch_multi_envelopes.r
envs.group <- loadenv('envs499_K_3.rds')
g <- plot.envs.grouped(envs.group$Thomas, ncol=4)
ggsave('report_out/07_envs.bar.thomas.pdf', plot=g, width=5.5, height=3)
g <- plot.envs.grouped(envs.group$MatClust, ncol=4)
ggsave('report_out/08_envs.bar.matclust.pdf', plot=g, width=5.5, height=3)

# 3.6 - Branching points analysis ==========
# A nice plot of all patterns
g <- grid.arrange(grobs=Map(pppplot, data$ppp, data.branching$ppp, names(data$ppp)))
ggsave('report_out/09_patterns.branch.pdf', plot=g, width=5.5, height=5)

# Under the model hypotheses, the parent points (which we take to be the
# branching points in this case) should be CSR
envs.csr.branching <- grouped(csrenvs, data.branching, nsim=999)
g <- plot.envs.grouped(envs.csr.branching, ncol=4)
ggsave('report_out/10_csr.branching.pdf', plot=g, width=5.5, height=3)

# Compare npoints with fitted params times area
df.devis <- data.frame(
  id=1:length(data$ppp), truth=sapply(data.branching$ppp, npoints),
  thomas=fit.each.thomas['kappa',]*areas, matclust=fit.each.matclust['kappa',]*areas
)
g <- ggplot(df.devis, aes(x=reorder(id, -rowMeans(cbind(thomas, matclust)-truth)))) +
  geom_point(aes(y=truth, color='actual')) +
  geom_point(aes(y=thomas, color='Thomas'), shape=15) +
  geom_point(aes(y=matclust, color='MatClust'), shape=17) +
  labs(x='pattern #', y=expression(italic(kappa*abs(W)))) +
  scale_colour_manual(values=c(c("#444444"),c("#F7951F"),c("#0154A6"))) +
  theme(legend.position='bottom', legend.title=element_blank())
ggsave('report_out/11_deviations.pdf', plot=g, width=5, height=3)
