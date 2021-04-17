# This script runs A LOT of hypothesis tests in a non-interactive environment.
# (Useful for batch use and multiprocessing.)
args <- commandArgs(trailingOnly=TRUE)
nsim <- strtoi(args[1])
out <- paste('envelopes/envs', nsim, '_single.rds', sep='')

# Initialize and fit models ----
library(spatstat)
source('util.r')
data <- loaddata()
source('multiGET.r')

# This is important for reproducibility!
set.seed(012101)

# Problems with very low point counts leading to broken estimates.
# Behaviour changed in spatstat 2 experimentally:
spatstat.options(kppm.strict=F)

plan(multicore)
#handlers(handler_progress(':spin [:bar] :percent (:current/:total) in :elapsed(:tick_rate) ETA :eta'))
#handlers(global=T)

# Envelopes ----
doenv <- function(ppp, cluster) {
  cat('  Fitting model\n')
  fit <- kppm(ppp, cluster=cluster)
  multiGET.composite(list(ppp), fit, Gest, alpha=0.05, type='erl', nsim=nsim)
}
cat('Running MatClust...\n')
envs.matclust <- as.anylist(future_lapply(data$ppp, doenv, cluster='MatClust', future.seed=T, future.stdout=NA))
cat('Running Thomas...\n')
envs.thomas   <- as.anylist(future_lapply(data$ppp, doenv, cluster='Thomas', future.seed=T, future.stdout=NA))

saveRDS(list(Thomas=envs.thomas, MatClust=envs.matclust), file=out)
