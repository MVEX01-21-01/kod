library(spatstat)
library(GET)
library(gridExtra)
library(future.apply)
library(progressr)

# Helper function for rule determining how to recalculate significance for
# multiple testing. Currently, we use the Šidák correction, as suggested
# by Baddeley: https://stackoverflow.com/a/38176292/15077901
# Remember that conservatism contrary to what one might expect is unsafe, as we
# are looking to NOT reject the null hypothesis...
msignf <- function(ntests, alpha=0.05) 1 - (1-alpha)^(1/ntests)

# Helper functions to extract simulator from fit
match.r   <- function(fit) match.fun(paste('r', fit$internal$model, sep=''))
partial.r <- function(fit) {
  fun <- match.r(fit)
  function (fit, nsim, win) do.call(fun, as.list(c(fit$clustpar, fit$modelpar['mu'], nsim=nsim, win=list(win))))
}
# Helper function to extract statistic estimator from fit
match.est <- function(fit) switch(fit$info$fname, 'K'=Kest, 'g'=pcf)
# Helper function to simulate a fit
rFit <- function(fit, nsim, win) partial.r(fit)(fit, nsim, win)

# Helper function to obtain range from statistic
rrange <- function(obs) {
  argname <- fvnames(obs, '.x')
  
  r <- obs[[argname]]
  alims <- attr(obs, 'alim')
  list(r=r, rmin=alims[1], rmax=alims[2])
}

#' Generate multiple global envelopes for a set of patterns, under a composite
#' null hypothesis represented by a model fit, at significance level alpha.
#' The test statistic is determined by stat.
#' @param X Set of patterns to test agains
#' @param fit Fit representing null hypothesis, type minconfit
#' @param stat Function calculating statstic to test with from ppp
#' @param alpha Significance level for envelopes
#' @param type Type of global envelope, according to GET types
#' @param nsim Number of subtests n
#' @param nsim2 Number of simulations s-1
#' @param raw Whether to pass through spatstat envelope logic
#' @return A global envelope
multiGET.composite <- function(X, fit, stat, alpha=0.05, type='erl', nsim=NULL, nsim2=NULL, raw=F) {
  # significance required for each
  gamma <- msignf(length(X), alpha)
  if (is.null(nsim)) {
    nsim <- ceiling(1/gamma - 1)
  }
  if (is.null(nsim2)) {
    nsim2 <- nsim
  }
  
  p <- progressor(along=X)
  rfit <- partial.r(fit)
  future_lapply(X, function(ppp) {
    obs <- stat(ppp)
    range <- rrange(obs)
    # 1. estimates already obtained through fit
    # 2. simulate null curve replicates
    if (!raw) {
      enve2 <- envelope(ppp, stat, nsim=nsim2, simulate=rfit(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
    } else {
      valname <- fvnames(obs, '.y')
      sims2 <- rfit(fit, nsim2, ppp)
      stats2 <- sapply(sims2, function(i) stat(i, r=range$r)[[valname]])
      enve2 <- create_curve_set(list(r=range$r, obs=obs[[valname]], sim_m=stats2, theo=obs[['theo']]))
    }
    
    # 3. simulate and estimate null replicate parameters
    sims3 <- rfit(fit, nsim, ppp)
    fits3 <- Map(kppm, X=sims3, cluster=fit$internal$model, statistic=fit$info$fname)
    
    # 4. simulate composite curves
    if (!raw) {
      enve4 <- Map(function(ppp, fit) {
        envelope(ppp, stat, nsim=nsim2, simulate=rfit(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
      }, ppp=sims3, fit=fits3)
    } else {
      enve4 <- Map(function(ppp, fit) {
        obs <- stat(ppp)
        sims4 <- rfit(fit, nsim2, ppp)
        stats4 <- sapply(sims4, function(i) stat(i, r=range$r)[[valname]])
        enve4 <- create_curve_set(list(r=range$r, obs=obs[[valname]], sim_m=stats4, theo=obs[['theo']]))
      }, ppp=sims3, fit=fits3)
    }
    
    
    # 5-7. construct envelopes
    p()
    GET.composite(X=enve2, X.ls=enve4, r_min=range$rmin, r_max=range$rmax, type=type, alpha=gamma)
  }, future.seed=T)
}

# Retrieve the minimum p (closest to rejection) in each group
minp <- function(envs) min(sapply(envs, function(e) attr(e, "p")))

# Show all envelopes in all groups
multiGET.plot <- function(tests) {
  grid.arrange(grobs=mapply(function(group, name) {
    p <- 1 - (1 - minp(group))^length(group)
    arrangeGrob(grobs=mapply(function(env, who) plot(env) + labs(caption=who),
                             group, names(group), SIMPLIFY=F), top=paste(name, ": p =", p))
  }, tests, names(tests), SIMPLIFY=F), ncol=length(tests))
}

