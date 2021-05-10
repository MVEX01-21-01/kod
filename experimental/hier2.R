library(spatstat)
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
      #max(dists)
    })))
  })))
}

mat.scale.est2 <- function(children, branches) {
  unlist(lapply(1:length(children), function(i) {
    x <- as.matrix(children[[i]])
    branch <- branches[[i]]
    trees <- unique(x[, 3])
    trees <- trees[trees %in% branch$Tree]
    max_dists <- unlist(lapply(trees, function(j) {
      coords <- x[x[, 3] == j, , drop=F]
      parent <- branch[branch$Tree == j, ]
      pds <- pairdist(rbind(coords[, 1:2], parent[1:2]))
      dists <- pds[dim(pds)[1], ]
      #max(dists)
    }))
  }))
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
  sqrt(sum(dists^2)/(2*length(dists)))
}

thomas.scale.est.alt <- function(children, branches) {
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
  N <- length(dists)
  sigma_hat <- sqrt(sum(dists^2)/(2*length(dists)))
  nom <- 4^N * factorial(N) * factorial(N - 1) * sqrt(N)
  denom <- factorial(2 * N) * sqrt(pi)
  sigma <- sigma_hat * nom / denom
  sigma_hat
}




thomas.scale.est2 <- function(children, branches) {
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
  dists
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


sim.mat.parents <- function(mu, scale, parents, nsim, window, kappa) {
  if (class(parents) == 'ppp') {
    parlist <- data.frame(X=parents$x, Y=parents$y, Tree=1:length(parents$x))
  } else {
    parlist <- parents
  }
  win_dilated <- dilation(window, scale * 1.01)
  win_complement <- complement.owin(w=window, frame=win_dilated)
  poispar <- rpoispp(lambda=kappa, win=win_complement)
  if (poispar$n != 0) {
    max_tree <- max(parlist$Tree)
    poisparlist <- data.frame(X=poispar$x, Y=poispar$y, Tree=(max_tree + 1):(max_tree + poispar$n))
    parlist <- rbind(parlist, poisparlist)
  }
  npar <- dim(parlist)[1]
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

sim.thom.parents <- function(mu, scale, parents, nsim, window, kappa, scale.dilation=4) {
  if (class(parents) == 'ppp') {
    parlist <- data.frame(X=parents$x, Y=parents$y, Tree=1:length(parents$x))
  } else {
    parlist <- parents
  }
  win_dilated <- dilation(window, scale * scale.dilation)
  win_complement <- complement.owin(w=window, frame=win_dilated)
  poispar <- rpoispp(lambda=kappa, win=win_complement)
  if (poispar$n != 0) {
    max_tree <- max(parlist$Tree)
    poisparlist <- data.frame(X=poispar$x, Y=poispar$y, Tree=(max_tree + 1):(max_tree + poispar$n))
    parlist <- rbind(parlist, poisparlist)
  }
  npar <- dim(parlist)[1]
  res <- list()
  for (j in 1:npar) {
    parent <- parlist[j, , drop=F]
    n <- rpois(nsim, mu)
    centre = c(parent$X, parent$Y)
    disc <- lapply(n, function(x) {
      d <- gausdisc(n=x, sigma=scale, centre=centre)
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

thomas.scale.est(data_moderate, data_moderate_b)
thomas.scale.est.alt(data_moderate, data_moderate_b)

mat_mod <- mat.scale.est2(data_moderate, data_moderate_b)
mat_mod <- mat_mod[mat_mod > 0]
mat_norm <- mat.scale.est2(data_normal, data_normal_b)
mat_norm <- mat_norm[mat_norm > 0]
qq_mod <- quantile(runif(100, max = mat.scale.est(data_moderate, data_moderate_b)), probs = seq(0.01,0.99,0.01))
qq_norm <- quantile(runif(100, max = mat.scale.est(data_normal, data_normal_b)), probs = seq(0.01,0.99,0.01))

par(mfrow=c(2,2))
qqmat_mod <- quantile(mat_mod, probs = seq(0.01,0.99,0.01))
qqmat_norm <- quantile(mat_norm, probs = seq(0.01,0.99,0.01))

plot(qq_mod, ylab='', main='Matèrn MODERATE')
points(qqmat_mod, pch='*', col='blue')

plot(qq_norm, ylab='', main='Matèrn NORMAL')
points(qqmat_norm, pch='*', col='blue')

qq2_norm <- quantile(abs(rnorm(100, sd = thomas.scale.est(data_normal, data_normal_b))), probs = seq(0.01,0.99,0.01))
qq2_mod <- quantile(abs(rnorm(100, sd = thomas.scale.est(data_moderate, data_moderate_b))), probs = seq(0.01,0.99,0.01))

qqthom_mod <- quantile(thomas.scale.est2(data_moderate, data_moderate_b), probs = seq(0.01,0.99,0.01))
qqthom_norm <- quantile(thomas.scale.est2(data_normal, data_normal_b), probs = seq(0.01,0.99,0.01))

plot(qq2_mod, qqthom_mod, ylab='', main='Thomas MODERATE')
points(qqthom_mod, pch='*', col='blue')

plot(qq2_norm, qqthom_norm, ylab='', main='Thomas NORMAL')
points(qqthom_norm, pch='*', col='blue')
par(mfrow=c(1,1))
