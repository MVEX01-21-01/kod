## "Echelon" / direct model

# Parameter estimators ----

# Based on tree data, estimate R (Matérn cluster radius)
echelon.est.R <- function(X) {
  #return(2*mean(nndist(X)))
  P <- attr(X, 'parents')
  trees <- attr(X, 'parentid')
  if(is.null(P)) stop('Missing parent data from pattern')
  
  # the ML estimator is obviously the tightest possible radius.
  # but it seems that means perform much better
  D <- sapply(1:max(trees), function(t) {
    if (t <= npoints(P)) {
      # if we know THE parent point,
      #   we can just look at the distances
      x <- X$x[trees == t] - P$x[t]
      y <- X$y[trees == t] - P$y[t]
      mean(sqrt(x^2 + y^2))
    } else {
      # otherwise, guess by looking at distances between points
      d <- dist(matrix(c(X$x[trees == t], X$y[trees == t]), ncol=2))
      mean(d)
    }
  })
  mean(D, na.rm=T)
}

# TODO
echelon.est.sigma <- function(X) {}

# Based on tree data, estimate kappa (parent point intensity)
#  There are two methods, since the data varies:
#  1) use the number of tagged clusters
#  2) use the observed parent points
#  Under the assumption of CSR, parents, the actual locations are irrelevant,
#  just the number of clusters.
echelon.est.kappa <- function(X, by.parents=F) {
  if (!by.parents) {
    max(attr(X, 'parentid'))/area(X)
  } else {
    intensity(attr(X, 'parents'))
  }
}

# For mu, the natural estimate works just fine
echelon.est.mu <- function(X, kappa=NULL, by.parents=F) {
  if (is.null(kappa)) kappa <- echelon.est.kappa(X,by.parents=by.parents)
  intensity(X)/kappa
}

# High-level functions ----
matclust.estEchelon <- function(X, by.parents=F) {
  kappa <- echelon.est.kappa(X, by.parents=by.parents)
  R <- echelon.est.R(X)
  mu <- echelon.est.mu(X, kappa=kappa)
  c(kappa=kappa, mu=mu, scale=R)
}

thomas.estEchelon <- function(X) {
  # TODO
}

# Simulate from echelon fit (TODO Thomas)
rEchelon <- function(par, nsim=1, win=owin(c(0,1),c(0,1))) {
  # TODO do it manually????
  rMatClust(
    kappa=par['kappa'],
    mu=par['mu'],
    scale=par['scale'],
    nsim=nsim,
    win=as.owin(win),
    saveparents=T
  )
}

# Main ====
library(spatstat)
library(GET)
library(magrittr)
source('util.r')
source('multiGET.r')
source('plotting.r')
data <- loaddata.full()

# K
matclust.K <- function(r, kappa, R) {
  z <- r/(2*R)
  h <- ifelse(z > 1, 1, 2 + (1/pi) * ((8 * z^2 - 4) * acos(z) - 2 * asin(z) + 4 * z * sqrt((1 - z^2)^3) - 6 * z * sqrt(1 - z^2)))
  pi * r^2 + h/kappa
}

tf <- matclust.estEchelon(data$ppp[[10]])
tK <- Kest(data$ppp[[10]], correction='isotropic')
plot(tK, sqrt(iso/pi)-r~r)
lines(tK$r, sqrt(matclust.K(tK$r, tf['kappa'], tf['scale'])/pi) - tK$r, col='red')

# envelopes
echenv <- function(ppp, stat, nsim, alpha, type='erl') {
  cat('  Fitting model\n')
  fit <- matclust.estEchelon(ppp)
  
  cat(paste(' Preparing individual envelope, with N=', nsim, '\n'))
  rfit <- rEchelon
  nsim2 <- nsim
  gamma <- alpha
  
  #{
    tic <- proc.time()
    obs <- stat(ppp)
    range <- rrange(obs)
    # 1. estimates already obtained through fit
    # 2. simulate null curve replicates
    enve2 <- envelope(ppp, stat, nsim=nsim2, simulate=rfit(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
    
    # 3. simulate and estimate null replicate parameters
    sims3 <- rfit(fit, nsim, ppp)
    fits3 <- Map(function(i) {
      # here bad patterns crash. replace those!
      repeat {
        nfit <- try(matclust.estEchelon(sims3[[i]], by.parents=F))
        if (class(nfit) != 'try-error' && nfit['scale'] > 0) return(nfit)
        cat(paste('  ...retrying', i, '(n=', npoints(sims3[[i]]), ')\n'))
        sims3[[i]] <- rfit(fit, 1, ppp)[[1]] # replace...
      }
    }, i=1:length(sims3))
    
    # 4. simulate composite curves
    enve4 <- Map(function(ppp, fit) {
      envelope(ppp, stat, nsim=nsim2, simulate=rfit(fit, nsim2, ppp), r=range$r, savefuns=T, verbose=F)
    }, ppp=sims3, fit=fits3)
    
    # 5-7. construct envelopes
    delta <- proc.time() - tic
    cat(paste('  Done in', delta[3] / 60, 'min, at', Sys.time(), '\n'))
    GET.composite(X=enve2, X.ls=enve4, type=type)#, alpha=gamma)
  #}
}

#plan(multicore, workers=4)
#set.seed(012101)
envsE.matclust <- anylapply(data$ppp, function(X) envelope(X, Lest, nsim=499, simulate=rEchelon(matclust.estEchelon(X), 499, X)), global=T)
names(envsE.matclust) <- names(data$ppp)
plot.envs.single(envsE.matclust)

cat('Running MatClust (echelon)...\n')
envsE.matclust <- as.anylist(future_Map(echenv, ppp=data$ppp, stat=list(Lest), nsim=199, alpha=0.05, future.seed=T, future.stdout=NA))
names(envsE.matclust) <- names(data$ppp)
plot.envs.single(envsE.matclust)

#saveRDS(list(MatClust=envsE.matclust), file='ecenvs.rds')

