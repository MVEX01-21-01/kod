library(spatstat)

data_moderate <- readRDS('DATA_ENFS/CALF_MODERATE_df')
data_normal   <- readRDS('DATA_ENFS/CALF_NORMALS_df')
data_moderate_b <- readRDS('DATA_ENFS/CALF_MODERATE_BRANCHING_df')
data_normal_b   <- readRDS('DATA_ENFS/CALF_NORMALS_BRANCH_df')
load("fits.RData")

data <-   hyperframe(
  g = factor(rep.int(c('MODERATE', 'BRANCH'), c(length(data_moderate), length(data_moderate_b)))),
  ppp = c(data_moderate, data_moderate_b)
)

r11seq <- seq(20,60,by=5)
r12seq <- seq(5,15,by=5)
r22seq <- seq(20,60,by=5)
s <- expand.grid(r11=r11seq, r12=r12seq, r22=r22seq)

pfn <- function(r11, r12, r22) {
  HierStrauss(radii=matrix(c(r11,r12,r12,r22),2,2), archy=c(2,1))
}

pfn2 <- function(r11, r12, r22) {
  HierStraussHard(iradii=matrix(c(r11,r12,r12,r22),2,2), archy=c(2,1))
}
sd <- as.solist(lapply(1:7, function(i) {superimpose(nerv=data_moderate[[i]], branch=data_moderate_b[[i]])}))
  
H <- hyperframe(Y = sd)

d_1 <- superimpose(nerv=data_moderate[[1]],branch=data_moderate_b[[1]])

fit.hierstrausshard <- function(data=data, explore=F) {
  
  r11seq <- seq(20,60,by=5)
  r12seq <- seq(5,15,by=5)
  r22seq <- seq(20,60,by=5)
  s <- expand.grid(r11=r11seq, r12=r12seq, r22=r22seq)
  
  prf <- tryCatch({
    profilepl(s, pfn, data ~ marks, correction="translate", fast=TRUE, verbose=T)
  }, error = function(e) {
      print(e)
    })
  if (explore) {
    popt <- as.numeric(prf$param[prf$iopt, , drop = FALSE])
    r11seq <- seq(popt[1] - 5,popt[1] + 5,by=1)
    r12seq <- seq(popt[2] - 5,popt[2] + 5,by=1)
    r22seq <- seq(popt[3] - 5,popt[3] + 5,by=1)
    s <- expand.grid(r11=r11seq, r12=r12seq, r22=r22seq)
    prf <- profilepl(s, pfn, data ~ marks, correction="translate", fast=TRUE, verbose=T)
  }
  popt <- as.numeric(prf$param[prf$iopt, , drop = FALSE])
  imat = matrix(c(popt[1], popt[2], popt[2], popt[3]), 2, 2)
  ppm(data ~ marks, HierStraussHard(iradii=imat, archy=c(2,1)))
}

fit.hierstrausshard.mppm <- function(data=data, explore=F) {
  
  r11seq <- seq(20,60,by=5)
  r12seq <- seq(5,15,by=5)
  r22seq <- seq(20,60,by=5)
  s <- expand.grid(r11=r11seq, r12=r12seq, r22=r22seq)
  
  prf <- tryCatch({
    profilepl.mppm(s, pfn, Y ~ marks, data=data, correction="translate", fast=TRUE, verbose=T)
  }, error = function(e) {
    print(e)
  })
  if (explore) {
    popt <- as.numeric(prf$param[prf$iopt, , drop = FALSE])
    r11seq <- seq(popt[1] - 5,popt[1] + 5,by=1)
    r12seq <- seq(popt[2] - 5,popt[2] + 5,by=1)
    r22seq <- seq(popt[3] - 5,popt[3] + 5,by=1)
    s <- expand.grid(r11=r11seq, r12=r12seq, r22=r22seq)
    prf <- profilepl.mppm(s, pfn, Y ~ marks, data=data, correction="translate", fast=TRUE, verbose=T)
  }
  popt <- as.numeric(prf$param[prf$iopt, , drop = FALSE])
  imat = matrix(c(popt[1], popt[2], popt[2], popt[3]), 2, 2)
  mppm(Y ~ marks, data=data, HierStraussHard(iradii=imat, archy=c(2,1)))
}

fit.hierstrausshard.mppm(H)

fits <- lapply(1:7, function(i) { fit.hierstrausshard(data=superimposed_data[[i]], T) })

avg_rmat <- Reduce('+', lapply(1:7, function(i) { as.matrix(fits[[i]]$interaction$par$iradii)}))/length(fits)
avg_hmat <- Reduce('+', lapply(1:7, function(i) { as.matrix(fits[[i]]$interaction$par$hradii)}))/length(fits)
avg_beta <- Reduce('+', lapply(1:7, function(i) {  s <- summary.ppm(fits[[i]], quick="no variances")
                                        s$trend$value}))/length(fits)
avg_gamma <- matrix(Reduce('+', lapply(1:7, function(i) {
  s <- summary.ppm(fits[[i]], quick="no variances")$interaction$printable
  s <- mapply(s, FUN=as.numeric)
  s[3] = s[2]
  s}))/length(fits), nrow=2, ncol=2)

f <- mppm(Y ~ marks, data=H, HierStraussHard(iradii=avg_rmat, archy=c(2,1)))
res <- residuals(f, type="Pearson")
smor <- with(hyperframe(res=res), Smooth(res, sigma=4))

M <- rmhmodel(cif="straushm", par=list(beta=avg_beta,gamma=avg_gamma,iradii=avg_rmat,hradii=avg_hmat),types=c("nerv","branch"), w=square(1000))
X <- rmh(M)
