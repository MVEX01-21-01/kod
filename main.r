### MVEX01-21-01
### This file should render the results of the report exactly,
### for reproducibility.

# Initialization ==========
if (!dir.exists('DATA_ENFS')) {
  stop('Can\'t find pattern subdir DATA_ENFS. Check your working directory.')
}

library(ggplot2)
theme_set(theme_minimal())
library(spatstat)

source('util.r')
source('repcluster.r')
source('multiGET.r')
source('plotting.r')
data <- loaddata()

# Clean output directory
unlink('report_out', recursive=T)
dir.create('report_out')

# Helper to unnest single envelopes
unlist1 <- function(l) lapply(l, function(i) i[[1]])


# CSR test on patterns ==========
envs.csr <- grouped(csrenvs, data, nsim=499)
g <- plot.envs.grouped(envs.csr)
ggsave('report_out/csr.pdf', plot=g, width=6, height=5)

# Individual parameter fits ==========
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
g.indparams <- ggplot(data.params, aes(group, value)) +
  geom_boxplot() + facet_wrap(~ model + param, scales='free')
ggsave('report_out/ind.box.pdf', plot=g.indparams, width=6, height=6)

# Individual envelopes ==========
# Load pre-generated envelope, from batch_single_envelopes.r
# TODO fix data format?
envs.ind <- loadenv('envs199_single.rds')
g <- plot.envs.single(unlist1(envs.ind$Thomas))
ggsave('report_out/envs.ind.Thomas.pdf', plot=g, width=5, height=5)
g <- plot.envs.single(unlist1(envs.ind$MatClust))
ggsave('report_out/envs.ind.MatClust.pdf', plot=g, width=5, height=5)

# Group parameter fits ==========
## Show the variance in individual K, and the fitted curves
K <- grouped(Kbar, data, correction='iso')
fit.thomas.K <- repcluster.estK(data, 'Thomas')
fit.matclust.K <- repcluster.estK(data, 'MatClust')
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
ggsave('report_out/bar.L.fit.pdf', plot=g, width=6, height=3.5)

## Overlay the group results on the boxplots
df.fit <- rbind(
  longparams(rbind(sapply(fit.thomas.K, function(f) f$modelpar)), 'Thomas', c('MODERATE','NORMAL')),
  longparams(rbind(sapply(fit.matclust.K, function(f) f$modelpar)), 'MatClust', c('MODERATE','NORMAL'))
)
df.fit[df.fit=='sigma' | df.fit=='R'] <- 'scale'
g <- g.indparams + geom_point(data=df.fit, shape=3, color='red')
ggsave('report_out/bar.box.pdf', plot=g, width=6, height=6)

# Group envelopes ==========
# Load pre-generated envelope, from batch_multi_envelopes.r
envs.group <- loadenv('envs499_K_2.rds')
g <- plot.envs.grouped(envs.group$Thomas)
ggsave('report_out/envs.bar.thomas.pdf', plot=g, width=5, height=5)
g <- plot.envs.grouped(envs.group$MatClust)
ggsave('report_out/envs.bar.matclust.pdf', plot=g, width=5, height=5)

# 3.6 - Branching points analysis ==========
# TODO for now, see branch_points.r

