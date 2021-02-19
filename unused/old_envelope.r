# Envelope test ----
# Suggestions on how to perform multiple testing:
# Baddeley: https://stackoverflow.com/a/38176292/15077901
# https://research.csiro.au/software/wp-content/uploads/sites/6/2015/02/Rspatialcourse_CMIS_PDF-Standard.pdf
envelopes <- function(X, groups, fits, statistic, alpha=0.05, global=F, ...) {
  as.anylist(mapply(function(X, fit) {
    # significance required for each
    gamma <- 1 - (1-alpha)^(1/length(X))
    nsim  <- ceiling(1/gamma - 1)
    
    # simulate from fit
    rFit <- match.fun(paste('r', fit$internal$model, sep=''))
    
    # construct envelopes
    anylapply(X, function(ppp) {
      sims <- do.call(rFit,
                      # Note that a "global" envelope requires double the simulations
                      as.list(c(fit$clustpar, fit$modelpar['mu'], nsim=(1+global)*nsim, win=list(ppp))))
      envelope(ppp, statistic, nsim=nsim, simulate=sims, global=global, ...)
    })
  }, X=split(X, groups), fit=fits))
}

# Looking at the minimum p * no of tests. This was a suggestion to try to
# evaluate the significance levels, but it's not even clear this value is meaningful
tests.thomas <- lapply(envs.thomas, function(X) lapply(X, mad.test))
sapply(tests.thomas, function(X) length(X)*min(sapply(X, function(t) t$p.value)))
tests.matclust <- lapply(envs.matclust, function(X) lapply(X, mad.test))
sapply(tests.matclust, function(X) length(X)*min(sapply(X, function(t) t$p.value)))


###
# Helper function to simulate from fit into curve set, based on window
simulatecurves <- function(obs, fit, nsim, win, stat) {
  argname <- fvnames(obs, '.x')
  valname <- fvnames(obs, '.y')
  #r <- fit$objargs$rvals
  sims <- rFit(fit, nsim, win)
  create_curve_set(list(obs=obs[[valname]], r=obs[[argname]],
                        sim_m=t(sapply(sims, function(i) stat(i)[[valname]]))))
}