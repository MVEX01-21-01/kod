# This script runs A LOT of hypothesis tests in a non-interactive environment.
# (Useful for batch use and multiprocessing.)

# Initialize and fit models ----
library(spatstat)
source('load_data.r')
source('repcluster.r')
#fit.thomas   <- repcluster.estK(data, 'Thomas')
#fit.matclust <- repcluster.estK(data, 'MatClust')
fit.thomas   <- repcluster.estpcf(data, 'Thomas')
fit.matclust <- repcluster.estpcf(data, 'MatClust')
out <- 'envs499_pcf.rds' #'envs499.rds'
source('multiGET.r')

plan(multicore)
handlers(handler_progress(':spin [:bar] :percent (:current/:total) in :elapsed(:tick_rate) ETA :eta'))
handlers(global=T)

# Envelopes ----
envs.thomas   <- grouped(multiGET.composite, data, fit.thomas, c(Gest), alpha=0.05, type='erl', nsim=499)
envs.matclust <- grouped(multiGET.composite, data, fit.matclust, c(Gest), alpha=0.05, type='erl', nsim=499)

saveRDS(list(Thomas=envs.thomas, MatClust=envs.matclust), file=out)
