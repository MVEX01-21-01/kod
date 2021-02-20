# This script runs A LOT of hypothesis tests in a non-interactive environment.
# (Useful for batch use and multiprocessing.)

# Initialize and fit models ----
library(spatstat)
source('load_data.r')
source('repcluster.r')
fit.thomas   <- repcluster.estK(data, 'Thomas')
fit.matclust <- repcluster.estK(data, 'MatClust')
source('multiGET.r')

# Envelopes ----
plan(multicore)
handlers(list(handler_progress(':spin [:bar] :percent (:current/:total) in :elapsed(:tick_rate) ETA :eta')))
with_progress({
  envs.thomas   <- grouped(multiGET.composite, data, fit.thomas, c(Gest), alpha=0.05, type='erl', nsim=299)
  envs.matclust <- grouped(multiGET.composite, data, fit.matclust, c(Gest), alpha=0.05, type='erl', nsim=299)
})

saveRDS(list(Thomas=envs.thomas, MatClust=envs.matclust), file = 'envs.rds')
