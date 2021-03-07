# Explore branching points
# Note: check working directory
library(ggplot2)
theme_set(theme_minimal())
library(spatstat)
source('util.r')
data <- loaddata()
data.branching <- loaddata.branching()

# Visualization ================================================================
grid.arrange(grobs=lapply(1:nrow(data), function(i) {
  # extract
  child.x <- data$ppp[[i]]$x
  child.y <- data$ppp[[i]]$y
  parent.x <- data.branching$ppp[[i]]$x
  parent.y <- data.branching$ppp[[i]]$y
  xlim <- data$ppp[[i]]$window$xrange
  ylim <- data$ppp[[i]]$window$yrange
  
  ggplot() + geom_point(aes(x=child.x, y=child.y), shape=1) +
    geom_point(aes(x=parent.x, y=parent.y), size=3, shape=4) +
    labs(title=NULL, tag=i, x=NULL, y=NULL) + coord_fixed(xlim=xlim, ylim=ylim, expand=F, clip='off') +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.border = element_rect(color='gray', fill=NA, size=0.5))
}))

# CSR ==========================================================================
# Under the model hypotheses, the parent points (which we take to be the
# branching points in this case) should be CSR
csr.envs <- grouped(function(Xs) {
  lapply(Xs, function(X) {
    env <- envelope(X, nsim=999, savefuns=T, fun=Lest,
                    simulate=expression(runifpoint(ex=X)), transform=expression(.-r),
                    correction='trans')
    global_envelope_test(env,type='erl')
  })
}, data.branching)

groupplot(csr.envs)

# Compare direct kappa estimates with fitted params ============================
deviations <- function(Xs, Refs, cluster) {
  ests <- params.each(Xs, cluster)['kappa',]
  truth <- sapply(Refs, intensity)
  data.frame(truth, ests, deltas=(ests-truth))
}
deviations.gr <- function(Xs, Refs, cluster) {
  d <- deviations(Xs, Refs, cluster)
  ggplot(d, aes(x=truth, y=deltas)) + geom_point() + geom_abline(slope=1,intercept=0)
}

deviations.gr(data$ppp, data.branching$ppp, 'Thomas')
