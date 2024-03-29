# Data loading
loaddata <- function() {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL')
  finishdata(data_moderate, data_normal)
}
loaddata.branching <- function() {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL_BRANCHING')
  finishdata(data_moderate, data_normal)
}
finishdata <- function(mod, nor) {
  hyperframe(
    g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(mod), length(nor)))),
    ppp = anylapply(c(mod, nor),addunit)
  )
}
loaddata.full <- function() {
  Xs <- loaddata()
  Ps <- loaddata.branching()
  moderate_df <- readRDS('DATA_ENFS/CALF_MODERATE_df')
  normal_df <- readRDS('DATA_ENFS/CALF_NORMALS_df')
  moderate_b_df <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING_df')
  normal_b_df <- readRDS('DATA_ENFS/CALF_NORMALS_BRANCH_df')
  Xt <- lapply(c(moderate_df, normal_df), function(df) df$Tree)
  Pt <- lapply(c(moderate_b_df, normal_b_df), function(df) df$Tree)
  
  # TODO: why arent the branching points complete data-wise?
  # We will have to use the parent patterns, but renumber the trees
  # so the parentids are correct.
  # Number the parent trees according to index, and roll over
  # the unknown trees to higher numbers
  Ys <- anylapply(1:length(Xs$ppp), function(i) {
    xtrees <- Xt[[i]]
    ptrees <- Pt[[i]]  # parent trees
    utrees <- setdiff(xtrees, ptrees)  # unknown trees
    
    remap <- sapply(1:max(xtrees, ptrees), function(i) {
      v <- match(i, ptrees)
      if (is.na(v)) v <- match(i, utrees) + length(ptrees)
      v
    })
    
    X <- Xs$ppp[[i]]
    attr(X, 'parents') <- Ps$ppp[[i]]
    attr(X, 'parentid') <- sapply(xtrees, function(t) remap[t])
    X
  })
  
  hyperframe(g = Xs$g, ppp = Ys)
}
loadenv <- function(env) readRDS(paste('envelopes/', env, sep=''))

# Helper function rescaling to proper unit
addunit <- function(ppp) rescale(ppp,1,'µm')

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

# Helper to unnest single envelopes
unlist1 <- function(l) lapply(l, function(i) i[[1]])
