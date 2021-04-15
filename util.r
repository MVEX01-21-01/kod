# Data loading
loaddata <- function(tree=F) {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL')
  h <- hyperframe(
    g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
    ppp = anylapply(c(data_moderate, data_normal),addunit)
  )
  if (tree) {
    data_moderate_df <- readRDS('DATA_ENFS/CALF_MODERATE_df')
    data_normal_df <- readRDS('DATA_ENFS/CALF_NORMALS_df')
    h$ppp <- Map(augmentTree, ppp=h$ppp, ext=c(data_moderate_df,data_normal_df))
  }
  h
}
loaddata.branching <- function(tree=F) {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL_BRANCHING')
  h <- hyperframe(
    g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
    ppp = anylapply(c(data_moderate, data_normal),addunit)
  )
  if (tree) {
    data_moderate_df <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING_df')
    data_normal_df <- readRDS('DATA_ENFS/CALF_NORMALS_BRANCH_df')
    h$ppp <- Map(augmentTree, ppp=h$ppp, ext=c(data_moderate_df,data_normal_df))
  }
  h
}
loadenv <- function(env) readRDS(paste('envelopes/', env, sep=''))

# Helper function rescaling to proper unit
addunit <- function(ppp) rescale(ppp,1,'µm')
# Helper function to augment ppps with Tree tags
augmentTree <- function(ppp, ext) setmarks(ppp, ext$Tree)

# Helper function splitting over groups.
# This enables us to write functions in terms of a single set of patterns,
# and then using this function compose to groups
#grouped <- function(f, data, ...) anylapply(split(data$ppp, data$g), f, ...)
grouped <- function(f, data, ..., SIMPLIFY=F) {
  as.anylist(mapply(f, split(data$ppp, data$g), ..., SIMPLIFY=SIMPLIFY))
}

# Helper function obtaining fit parameters for a list of patterns
params.each <- function(X, cluster) {
  sapply(Map(kppm, X=X, cluster=cluster), function(fit) c(fit$clustpar, mu=fit$mu))
}

# Helper function performing CSR tests using L(r) and
# GET logic on a list of ppps.
csrenvs <- function(Xs, nsim=999) {
  lapply(Xs, function(X) {
    env <- envelope(X, nsim=nsim, savefuns=T, fun=Lest,
                    simulate=expression(runifpoint(ex=X)), transform=expression(.-r),
                    correction='trans')
    global_envelope_test(env,type='erl')
  })
}
