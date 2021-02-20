# Helper function splitting over groups.
# This enables us to write functions in terms of a single set of patterns,
# and then using this function compose to groups
#grouped <- function(f, data, ...) anylapply(split(data$ppp, data$g), f, ...)
grouped <- function(f, data, ..., SIMPLIFY=F) {
  as.anylist(mapply(f, split(data$ppp, data$g), ..., SIMPLIFY=SIMPLIFY))
}

# Load data
data_moderate <- readRDS('ENDPOINT_DATA/CALF_MODERATE')
data_normal   <- readRDS('ENDPOINT_DATA/CALF_NORMAL')
data <- hyperframe(
  g = factor(rep.int(c('MODERATE', 'NORMAL'), c(length(data_moderate), length(data_normal)))),
  ppp = c(data_moderate, data_normal)
)
