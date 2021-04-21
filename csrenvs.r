# csr envelopes in report 
library(spatstat)
library(GET)
source('util.r')

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

data <- loaddata()
data.branching <- loaddata.branching()

# Main
# This is important for reproducibility!
set.seed(012101)
nsim <- 999

envs.csr <- csrenvs(data$ppp, nsim=nsim)
saveRDS(envs.csr, file='envelopes/csrenv_endpoints.rds')

envs.csr.branching <- grouped(csrenvs, data.branching, nsim=nsim)
saveRDS(envs.csr.branching, file='envelopes/csrenv_branching.rds')