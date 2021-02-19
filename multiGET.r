library(GET)
library(gridExtra)

# Helper function for rule determining how to recalculate significance for
# multiple testing. Currently, we use the Šidák correction, as suggested
# by Baddeley: https://stackoverflow.com/a/38176292/15077901
# Remember that conservatism contrary to what one might expect is unsafe, as we
# are looking to NOT reject the null hypothesis...
msignf <- function(ntests, alpha=0.05) 1 - (1-alpha)^(1/ntests)

# Helper function to extract simulator from fit
match.r   <- function(fit) match.fun(paste('r', fit$internal$model, sep=''))
# Helper function to extract statistic estimator from fit
match.est <- function(fit) switch(fit$info$fname, 'K'=Kest, 'g'=pcf)
# Helper function to simulate a fit
rFit <- function(fit, nsim, win, fun=NULL) {
  do.call(ifelse(is.null(fun), match.r(fit), fun),
          as.list(c(fit$clustpar, fit$modelpar['mu'], nsim=nsim, win=list(win))))
}

# Helper function to obtain range from statistic
rrange <- function(obs) {
  argname <- fvnames(obs, '.x')
  #valname <- fvnames(obs, '.y')
  r <- obs[[argname]]
  alims <- attr(obs, 'alim')
  list(r=r, rmin=alims[1], rmax=alims[2])
}

# Generate multiple global envelopes for a set of patterns, under a composite
# null hypothesis represented by a model fit, at significance level alpha.
# The test statistic is determined by stat.
multiGET.composite <- function(X, fit, stat, alpha=0.05, type='area') {
  # significance required for each
  gamma <- msignf(length(X), alpha)
  nsim  <- ceiling(1/gamma - 1)
  nsim2 <- nsim # TODO?
  
  anylapply(X, function(ppp) {
    fun <- match.r(fit)
    range <- rrange(stat(ppp))
    
    # 1. estimates already obtained through fit
    # 2. simulate null curve replicates
    enve2 <- envelope(ppp, stat, nsim=nsim2, simulate=rFit(fit, nsim2, ppp, fun=fun), r=range$r, savefuns=T)
    
    # 3. simulate and estimate null replicate parameters
    sims3 <- rFit(fit, nsim, ppp, fun=fun)
    fits3 <- mapply(kppm, X=sims3, cluster=fit$internal$model, statistic=fit$info$fname, SIMPLIFY=F)
    
    # 4. simulate composite curves
    enve4 <- mapply(function(ppp, fit) {
      envelope(ppp, stat, nsim=nsim2, simulate=rFit(fit, nsim2, ppp, fun=fun), r=range$r, savefuns=T)
    }, ppp=sims3, fit=fits3, SIMPLIFY=F)
    
    # 5-7. construct envelopes
    GET.composite(X=enve2, X.ls=enve4, r_min=range$rmin, r_max=range$rmax, type=type, alpha=gamma)
  })
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

