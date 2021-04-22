library(spatstat)
setwd("~/Kandidatarbete/kod")
source("util.r")
data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE_df')
data_normal   <- readRDS('DATA_ENFS/CALF_NORMALS_df')
data_moderate_b <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING_df')
data_normal_b   <- readRDS('DATA_ENFS/CALF_NORMALS_BRANCH_df')
data_moderate_notree <- readRDS('DATA_ENFS/CALF_MODERATE')
data_normal_notree <- readRDS('DATA_ENFS/CALF_NORMAL')

mat.scale.est <- function(children, branches) {
  max(unlist(lapply(1:length(children), function(i) {
    x <- as.matrix(children[[i]])
    branch <- branches[[i]]
    trees <- unique(x[, 3])
    trees <- trees[trees %in% branch$Tree]
    max_dists <- max(unlist(lapply(trees, function(j) {
      coords <- x[x[, 3] == j, , drop=F]
      parent <- branch[branch$Tree == j, ]
      pds <- pairdist(rbind(coords[, 1:2], parent[1:2]))
      dists <- pds[dim(pds)[1], ]
      max(dists)
    })))
  })))
}

thomas.scale.est <- function(children, branches) {
  dists <- unlist(lapply(1:length(children), function(i) {
    x <- as.matrix(children[[i]])
    branch <- branches[[i]]
    trees <- unique(x[, 3])
    trees <- trees[trees %in% branch$Tree]
    dists <- lapply(trees, function(j) {
      coords <- x[x[, 3] == j, 1:2, drop=F]
      parent <- branch[branch$Tree == j, 1:2]
      pds <- pairdist(rbind(coords, parent))
      pds[dim(pds)[1], 1:dim(pds)[2] - 1]
    })
  }))
  sqrt(sum(dists^2)/(length(dists) - 1))
}

mu.est <- function(children, branches) {
  bc <- colSums(matrix(unlist(lapply(1:length(children), function(i) {
    child <- children[[i]]
    branch <- branches[[i]]
    branches <- length(branch$Tree)
    children <- length(child$Tree[child$Tree %in% branch$Tree])
    c(branches, children)
  })), ncol = 2, byrow = T))
  bc[2] / bc[1]
}


sim.mat.parents <- function(mu, scale, parents, nsim, window) {
  if (class(parents) == 'ppp') {
    npar <- parents$n
    parlist <- data.frame(X=parents$x, Y=parents$y, Tree=1:length(parents$x))
  } else {
    npar <- dim(parents)[1]
    parlist <- parents
  }
  if (is.null(npar)) {
    stop('Invalid parents')
  }
  res <- list()
  for (j in 1:npar) {
    parent <- parlist[j, , drop=F]
    n <- rpois(nsim, mu)
    centre = c(parent$X, parent$Y)
    disc <- lapply(n, function(x) {
      d <- runifdisc(n=x, radius=scale, centre=centre)
      #Create new ppp to remove invalid points
      ppp <- ppp(d$x, d$y, window = window)
      attr(ppp, 'parentid') <- rep(parent$Tree, ppp$n)
      ppp
    })
    res[[j]] <- disc
  }
  superimposed.res <- solist()
  for (i in 1:nsim) {
    dd <- solapply(1:npar, function(k) {res[[k]][[i]]})
    pid <- unlist(lapply(1:npar, function(i) {attr(dd[[i]], 'parentid')}))
    superimposed.res[[i]] <- superimpose(dd, W=window)
    attr(superimposed.res[[i]], 'parents') <- parents
    attr(superimposed.res[[i]], 'parentid') <- pid
  }
  if (nsim == 1) {
    superimposed.res[[1]]
  } else {
    superimposed.res
  }
}

sim.thom.parents <- function(mu, scale, parents, nsim, window) {
  if (class(parents) == 'ppp') {
    npar <- parents$n
    parlist <- data.frame(X=parents$x, Y=parents$y, Tree=1:length(parents$x))
  } else {
    npar <- dim(parents)[1]
    parlist <- parents
  }
  if (is.null(npar)) {
    stop('Invalid parents')
  }
  res <- list()
  sigma <- scale / sqrt(2)
  for (j in 1:npar) {
    parent <- parlist[j, , drop=F]
    n <- rpois(nsim, mu)
    centre = c(parent$X, parent$Y)
    disc <- lapply(n, function(x) {
      d <- gausdisc(n=x, sigma=sigma, centre=centre)
      ppp <- ppp(d[, 1], d[, 2], window=window)
      attr(ppp, 'parentid') <- rep(parent$Tree, ppp$n)
      ppp
      })
    res[[j]] <- disc
  }
  superimposed.res <- solist()
  for (i in 1:nsim) {
    dd <- solapply(1:npar, function(k) {res[[k]][[i]]})
    pid <- unlist(lapply(1:npar, function(i) {attr(dd[[i]], 'parentid')}))
    superimposed.res[[i]] <- superimpose(dd, W=window)
    attr(superimposed.res[[i]], 'parents') <- parents
    attr(superimposed.res[[i]], 'parentid') <- pid
  }
  if (nsim == 1) {
    superimposed.res[[1]]
  } else {
    superimposed.res
  }
}

gausdisc <- function(n, sigma, centre) {
  res <- matrix(rnorm(2 * n, mean=0, sd=sigma), ncol=2)
  sweep(res, 2, centre, "+")
}