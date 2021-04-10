# This script runs A LOT of hypothesis tests in a non-interactive environment.
# (Useful for batch use and multiprocessing.)

args <- commandArgs(trailingOnly=TRUE)
nsim <- strtoi(args[1])
out <- paste('envelopes/envs', nsim, '.rds', sep='')

# Initialize and fit models ----
library(spatstat)
source('util.r')
data <- loaddata()
source('repcluster.r')
fit.thomas   <- repcluster.estK(data, 'Thomas')
fit.matclust <- repcluster.estK(data, 'MatClust')
#fit.thomas   <- repcluster.estpcf(data, 'Thomas')
#fit.matclust <- repcluster.estpcf(data, 'MatClust')
source('multiGET.r')

plan(multicore)
#handlers(handler_progress(':spin [:bar] :percent (:current/:total) in :elapsed(:tick_rate) ETA :eta'))
#handlers(global=T)

# Envelopes ----
message('Running Thomas...')
envs.thomas   <- grouped(multiGET.composite, data, fit.thomas, c(Gest), alpha=0.05, type='erl', nsim=nsim)
message('Running MatClust...')
envs.matclust <- grouped(multiGET.composite, data, fit.matclust, c(Gest), alpha=0.05, type='erl', nsim=nsim)

saveRDS(list(Thomas=envs.thomas, MatClust=envs.matclust), file=out)
