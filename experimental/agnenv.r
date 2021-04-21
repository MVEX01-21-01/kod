library(spatstat)
library(GET)
source('multiGET.r')

#' Generate multiple global envelopes for a set of patterns, under a composite
#' null hypothesis represented by a model fit, at significance level alpha.
#' This formulation is model-agnostic compared to multiGET.composite, but
#' simplified in its functionality.
#' To test envelopes for something, you need a function to fit a model to an
#' observation, a function to simulate given a fit, and a function to compute
#' the statistic given an observation.
#' @param X List of patterns to test
#' @param statfun Function calculating statistic to test with from ppp, returns fv
#' @param fitfun Function calculating general fit object from ppp
#' @param simfun Function (fit, nsim, ppp) simulating ppp from general fit object
#' @param fit Starting fit representing null hypothesis, general. If NULL, use fitfun to compute for single pattern
#' @param alpha Significance level for collective test
#' @param type Type of global envelope, according to GET types
#' @param nsim Number of subtests N
#' @param nsim2 Number of simulations M-1
#' @return A global envelope
agnenv.composite <- function(X, statfun, fitfun, simfun, fit=NULL, alpha=0.05, type='erl', nsim=NULL, nsim2=NULL) {
  # significance required for each
  gamma <- msignf(length(X), alpha)
  if (is.null(nsim)) {
    nsim <- ceiling(1/gamma - 1)
  }
  if (is.null(nsim2)) {
    nsim2 <- nsim
  }
  
  if (is.null(fit)) {
    if (length(X) == 1) {
      cat('  Fitting model...\n')
      fit <- fitfun(X[[1]])
    } else
      stop('No null model given for replicates (fit is NULL)')
  } 
  
  cat(paste(' Preparing', length(X), 'envelope(s), with N=', nsim, '\n'))
  future_lapply(X, function(ppp) {
    tic <- proc.time()
    obs <- statfun(ppp)
    range <- rrange(obs)
    # 1. estimates already obtained through fit
    # 2. simulate null curve replicates
    enve2 <- envelope(ppp, statfun, nsim=nsim2, simulate=simfun(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
  
    # 3. simulate and estimate null replicate parameters
    sims3 <- simfun(fit, nsim, ppp)
    fits3 <- Map(function(i) {
      # here bad patterns crash. replace those!
      repeat {
        nfit <- try(fitfun(sims3[[i]]))
        if (class(nfit) != 'try-error') return(nfit)
        cat(paste('  ...retrying', i, '(n=', npoints(sims3[[i]]), ')\n'))
        sims3[[i]] <- rfit(fit, 1, ppp)[[1]] # replace...
      }
    }, i=1:length(sims3))
  
    # 4. simulate composite curves
    enve4 <- Map(function(ppp, fit) {
      envelope(ppp, statfun, nsim=nsim2, simulate=simfun(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
    }, ppp=sims3, fit=fits3)
  
    # 5-7. construct envelopes
    delta <- proc.time() - tic
    cat(paste('  Done in', delta[3] / 60, 'min, at', Sys.time(), '\n'))
    print(curvesfinite(enve2))
    print(sapply(enve4, curvesfinite))
    GET.composite(X=enve2, X.ls=enve4, r_min=range$rmin, r_max=range$rmax, type=type, alpha=gamma)
  }, future.seed=T, future.stdout=NA)
}

# Helper function to check if envelope finite
curvesfinite <- function(envelope) {
  if (!all(is.finite(envelope$obs))) return(F)
  all(sapply(attr(envelope, 'simfuns'), function(sim) {
    all(is.finite(sim))
  }))
}