setwd("~/Kandidatarbete/kod")
source('experimental/agnenv.r')
source('experimental/hier2.R')
source('util.r')
source('plotting.r')

# The new loaddata.full() augments with tree data internally.
data <- loaddata.full()

# Example for Thomas model
fitfun.thomas <- function(ppp) {
  # Translate into your data structures.
  # Sadly a bit ineffective, but using the spatstat internal support for
  # parent data allows us to contain all info in a single ppp.
  pppb <- attr(ppp, 'parents')
  trees <- attr(ppp, 'parentid')
  if(is.null(pppb)) stop('Missing parent data from pattern')
  children <- list(data.frame(X=ppp$x, Y=ppp$y, Tree=trees))
  branches <- list(data.frame(X=pppb$x, Y=pppb$y, Tree=1:length(pppb$x)))
  result <- list(
    # Note: you might want to design your own kappa estimator. Here it uses
    # the number of trees over the area (ML-estimate for Poisson I think).
    kappa=length(unique(sort(trees)))/area(ppp)
  )
  # Muting the output (otherwise internal assignment would be better)
  capture.output(result$mu <- mu.est(children, branches))
  capture.output(result$scale <- thomas.scale.est(children, branches))
  
  # Make sure to reject patterns that are problematic in some way.
  # The agnenv function will retry them. Otherwise, if you output a garbage fit,
  # simfun or envelope will crash later, which is unrecoverable.
  if (result$scale <= 0) stop('bad scale')
  result
}

fitfun.mat <- function(ppp) {
  # Translate into your data structures.
  # Sadly a bit ineffective, but using the spatstat internal support for
  # parent data allows us to contain all info in a single ppp.
  pppb <- attr(ppp, 'parents')
  trees <- attr(ppp, 'parentid')
  if(is.null(pppb)) stop('Missing parent data from pattern')
  children <- list(data.frame(X=ppp$x, Y=ppp$y, Tree=trees))
  branches <- list(data.frame(X=pppb$x, Y=pppb$y, Tree=1:length(pppb$x)))
  result <- list(
    # Note: you might want to design your own kappa estimator. Here it uses
    # the number of trees over the area (ML-estimate for Poisson I think).
    kappa=length(unique(sort(trees)))/area(ppp)
  )
  # Muting the output (otherwise internal assignment would be better)
  capture.output(result$mu <- mu.est(children, branches))
  capture.output(result$scale <- mat.scale.est(children, branches))
  
  # Make sure to reject patterns that are problematic in some way.
  # The agnenv function will retry them. Otherwise, if you output a garbage fit,
  # simfun or envelope will crash later, which is unrecoverable.
  if (result$scale <= 0) stop('bad scale')
  result
}

simfun <- function(fit, nsim, ppp) {
  # Just extract the parameters and use spatstat functions, no need to do it
  # all by yourself. Save the parents in the ppp so they can be re-used later!
  rThomas(fit$kappa, fit$scale, fit$mu, win=as.owin(ppp), nsim=nsim, saveparents=T)
}

simfun.thomas <- function(fit, nsim, ppp) {
  pppb <- attr(ppp, 'parents')
  sim.thom.parents(mu=fit$mu, scale=fit$scale, parents=pppb, win=as.owin(ppp), nsim=nsim, kappa=fit$kappa)
}

simfun.mat <- function(fit, nsim, ppp) {
  pppb <- attr(ppp, 'parents')
  sim.mat.parents(mu=fit$mu, scale=fit$scale, parents=pppb, win=as.owin(ppp), nsim=nsim, kappa=fit$kappa)
}

# Run some envelopes ====
plan(multicore)
envs.ind.hier2.thom´as <- future_lapply(data$ppp, function(X) {
  agnenv.composite(list(X), Lest, fitfun.thomas, simfun.thomas, nsim=19)
}, future.seed=T, future.stdout=NA)

envs.ind.hier2.mat <- future_lapply(data$ppp, function(X) {
  agnenv.composite(list(X), Lest, fitfun.mat, simfun.mat, nsim=19)
}, future.seed=T, future.stdout=NA)

envs.ind.hier2.mat <- sapply(1:length(envs.ind.hier2.mat), function(i) {envs.ind.hier2.mat[i][[1]]})
envs.ind.hier2.thomas <- sapply(1:length(envs.ind.hier2.thomas), function(i) {envs.ind.hier2.thomas[i][[1]]})

saveRDS(envs.ind.hier2.thomas, 'envelopes/hier2env.rds')
