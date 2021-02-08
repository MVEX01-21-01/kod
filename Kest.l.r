# Function to estimate K for listof ppp
Kest.l <- function(lX, ...) {
  # number of points square
  nsq <- sapply(lX, function(ppp) ppp$n^2)
  
  # minimum rmax, according to the rule (see Badderley)
  rmax <- min(sapply(lX, function(ppp) rmax.rule('K', Window(ppp), intensity(ppp))))
  
  # Extract all K estimates
  Kests <- lapply(lX, function(ppp) Kest(ppp, rmax=rmax, ...))
  # Extract a theoretical function computed by spatstat
  theo <- Kests[[1]][, 2, drop=T]
  # Extract a combined function computed through weighting (see Myllymäki, Panoutsopoulou, and Särkkä)
  bar <- rowSums(mapply(function(K,n2) K[, 3, drop=T]*n2, Kests, nsq))/sum(nsq)
  
  # Combine into the resulting object
  do.call(cbind, lapply(Kests, function(K) K[, -2])) %>%
    tweak.fv.entry('iso',new.tag='iso.0') %>%
    bind.fv(data.frame(bar=bar), 'bar(%s)[iso](r)', 'weighted mean of hat(%s)(r)', 'bar') %>%
    bind.fv(data.frame(theo=theo), '%s[pois](r)', 'theoretical Poisson %s(r)', 'bar')
}