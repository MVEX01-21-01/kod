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

test2 <- thomas.scale.est(data_moderate, data_moderate_b)


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

mat_scale_moderate <- mat.scale.est(data_moderate, data_moderate_b)
mat_scale_normal <- mat.scale.est(data_normal, data_normal_b)
mu_moderate <- mu.est(data_moderate, data_moderate_b)
mu_normal <- mu.est(data_normal, data_normal_b)
thom_scale_moderate <- thomas.scale.est(data_moderate, data_moderate_b)
thom_scale_normal <- thomas.scale.est(data_normal, data_normal_b)


sim.mat.parents <- function(mu, scale, parents, nsim, window) {
  npar <- dim(parents)[1]
  res <- list()
  for (j in 1:npar) {
    parent <- parents[j, , drop=F]
    n <- rpois(nsim, mu)
    centre = c(parent$X, parent$Y)
    disc.ppp <- lapply(n, function(x) {runifdisc(n=x, radius=scale, centre=centre)})
    res[[j]] <- disc.ppp
  }
  superimposed.res <- list()
  for (i in 1:nsim) {
    dd <- solapply(1:npar, function(k) {res[[k]][[i]]})
    superimposed.res[[i]] <- superimpose(dd, W=window)
  }
  superimposed.res
}

sim.thom.parents <- function(mu, scale, parents, nsim, window) {
  npar <- dim(parents)[1]
  res <- list()
  sigma <- scale / sqrt(2)
  for (j in 1:npar) {
    parent <- parents[j, , drop=F]
    n <- rpois(nsim, mu)
    centre = c(parent$X, parent$Y)
    disc <- lapply(n, function(x) {
      d <- gausdisc(n=x, sigma=sigma, centre=centre)
      ppp(d[, 1], d[, 2], window=window)
      })
    res[[j]] <- disc
  }
  superimposed.res <- list()
  for (i in 1:nsim) {
    dd <- solapply(1:npar, function(k) {res[[k]][[i]]})
    superimposed.res[[i]] <- superimpose(dd, W=window)
  }
  superimposed.res
}

gausdisc <- function(n, sigma, centre) {
  res <- matrix(rnorm(2 * n, mean=0, sd=sigma), ncol=2)
  sweep(res, 2, centre, "+")
}


test <- sim.mat.parents(mu_moderate, mat_scale_moderate, data_moderate_b[[1, drop=F]], 10, data_moderate_notree[[1]]$window)
test2 <- sim.thom.parents(mu_moderate, mat_scale_moderate, data_moderate_b[[1, drop=F]], 10, data_moderate_notree[[1]]$window)
plot.ppp(test2[[7]])
