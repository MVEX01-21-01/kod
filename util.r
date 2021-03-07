# Data loading
loaddata <- function() {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL')
  hyperframe(
    g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
    ppp = c(data_moderate, data_normal)
  )
}
loaddata.branching <- function() {
  data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING')
  data_normal   <- readRDS('DATA_ENFS/CALF_NORMAL_BRANCHING')
  hyperframe(
    g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
    ppp = c(data_moderate, data_normal)
  )
}

# Helper function splitting over groups.
# This enables us to write functions in terms of a single set of patterns,
# and then using this function compose to groups
#grouped <- function(f, data, ...) anylapply(split(data$ppp, data$g), f, ...)
grouped <- function(f, data, ..., SIMPLIFY=F) {
  as.anylist(mapply(f, split(data$ppp, data$g), ..., SIMPLIFY=SIMPLIFY))
}

# Helper function obtaining fit parameters for a list of patterns
params.each <- function(X, cluster) {
  sapply(Map(kppm, X=X, cluster=cluster), function(fit) fit$clustpar)
}

# Helper function plotting nested and factored lists
groupplot <- function(groups, topfun=NULL) {
  grid.arrange(grobs=Map(function(group, name) {
    arrangeGrob(grobs=Map(function(x, who) plot(x) + labs(tag=who),
                          group, names(group)), top=switch(is.null(topfun)+1,paste(name, topfun(group, x)),name))
  }, groups, names(groups)), ncol=length(groups))
}
