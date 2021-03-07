# This script runs A LOT of hypothesis tests in a non-interactive environment.
# (Useful for batch use and multiprocessing.)
args <- commandArgs(trailingOnly=TRUE)
nsim <- strtoi(args[1])
out <- paste('envs', nsim, '_single.rds', sep='')

# Initialize and fit models ----
library(spatstat)
source('util.r')
data <- loaddata()
source('multiGET.r')

plan(multicore)
handlers(handler_progress(':spin [:bar] :percent (:current/:total) in :elapsed(:tick_rate) ETA :eta'))
handlers(global=T)

# Envelopes ----
doenv <- function(ppp, cluster) {
  fit <- kppm(ppp, cluster=cluster)
  multiGET.composite(list(ppp), fit, Gest, alpha=0.05, type='erl', nsim=nsim)
}
envs.thomas   <- as.anylist(future_lapply(data$ppp, doenv, cluster='Thomas'))
envs.matclust <- as.anylist(future_lapply(data$ppp, doenv, cluster='MatClust'))

saveRDS(list(Thomas=envs.thomas, MatClust=envs.matclust), file=out)
