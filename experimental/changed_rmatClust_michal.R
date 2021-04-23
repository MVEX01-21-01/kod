#simulate endpoints as runifdisc around parent points
rPMatClust<-function (parent,kappa=1, scale, mu, win = owin(c(0, 1), c(0, 1)), nsim = 1, 
                      drop = TRUE, saveLambda = FALSE, expand = scale, ..., poisthresh = 1e-06, 
                      saveparents = TRUE) 
{
  win=parent$win
  if (missing(scale)) 
    scale <- list(...)$r
  kok <- is.numeric(kappa) || is.im(kappa)
  if (kok) {
    kappamax <- max(kappa)
  }
  else {
    kim <- as.im(kappa, W = win, ..., strict = TRUE)
    kra <- range(kim)
    kappamax <- kra[2] + 0.05 * diff(kra)
  }
  if (1/(pi * kappamax * scale^2) < poisthresh) {
    kapmu <- mu * (if (kok) 
      kappa
      else kim)
    result <- rpoispp(kapmu, win = win, nsim = nsim, drop = drop, 
                      warnwin = FALSE)
    return(result)
  }
  result <- simAroundParents(parent,kappa, scale, list(mu, function(x){runifdisc(x,radius=scale)}), 
                             win, radius = scale, nsim = nsim, drop = FALSE, saveparents = saveparents || 
                               saveLambda)
  if (saveLambda) {
    for (i in 1:nsim) {
      parents <- attr(result[[i]], "parents")
      Lambda <- clusterfield("MatClust", parents, scale = scale, 
                             mu = mu, ...)
      attr(result[[i]], "Lambda") <- Lambda[win, drop = FALSE]
    }
  }
  return(if (nsim == 1 && drop) result[[1]] else result)
}
#function used by rPMatClust
simAroundParents<-function(parentlist,kappa, expand, rcluster, win = owin(c(0, 1), c(0, 
                                                         1)), ..., lmax = NULL, nsim = 1, drop = TRUE, nonempty = TRUE, 
          saveparents = TRUE)
  {
  if (missing(expand) && !is.null(rmax <- list(...)$rmax)) 
    expand <- rmax
  if (is.function(rcluster)) 
    return(rPoissonCluster(kappa, expand, rcluster, win, 
                           ..., lmax = lmax, nsim = nsim, drop = drop, saveparents = saveparents))
  if (!(is.list(rcluster) && length(rcluster) == 2)) 
    stop("rcluster should be either a function, or a list of two elements")
  win <- as.owin(win)
  mu <- rcluster[[1]]
  rdisplace <- rcluster[[2]]
  if (is.numeric(mu)) {
    if (!(length(mu) == 1 && mu >= 0)) 
      stop("rcluster[[1]] should be a single nonnegative number")
    mumax <- mu
  }
  else if (is.im(mu) || is.function(mu)) {
    if (is.function(mu)) 
      mu <- as.im(mu, W = win, ..., strict = TRUE)
    mumax <- max(mu)
  }
  else stop("rcluster[[1]] should be a number, a function or a pixel image")
  if (!is.function(rdisplace)) 
    stop("rcluster[[2]] should be a function")
  frame <- boundingbox(win)
  dilated <- grow.rectangle(frame, expand)
  if (is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa))) 
    stop(paste("The window in which the image", sQuote("kappa"), 
               "is defined\n", "is not large enough to contain the dilation of the window", 
               sQuote("win")))
  if (nonempty) {
    if (is.function(kappa)) {
      kappa <- as.im(kappa, W = dilated, ..., strict = TRUE)
      lmax <- NULL
    }
    kappa <- kappa * (1 - exp(-mumax))
  }
  resultlist <- vector(mode = "list", length = nsim)
  
  for (i in 1:nsim) {
    parents <- parentlist
    np <-parents$n 
    if (np == 0) {
      result <- ppp(numeric(0), numeric(0), window = win)
      parentid <- integer(0)
    }
    else {
      if (!nonempty) {
        csize <- rpois(np, mumax)
      }
      else {
        csize <- qpois(runif(np, min = dpois(0, mumax)), 
                       mumax)
      }
      noff <- sum(csize)
      xparent <- parents$x
      yparent <- parents$y
      x0 <- rep.int(xparent, csize)
      y0 <- rep.int(yparent, csize)
      dd <- rdisplace(noff)
      mm <- if (is.ppp(dd)) 
        marks(dd)
      else NULL
      xy <- xy.coords(dd)
      dx <- xy$x
      dy <- xy$y
      if (!(length(dx) == noff)) 
        stop("rcluster returned the wrong number of points")
      xoff <- x0 + dx
      yoff <- y0 + dy
      parentid <- rep.int(1:np, csize)
      retain <- inside.owin(xoff, yoff, win)
      if (is.im(mu)) 
        retain[retain] <- inside.owin(xoff[retain], 
                                      yoff[retain], as.owin(mu))
      xoff <- xoff[retain]
      yoff <- yoff[retain]
      parentid <- parentid[retain]
      if (!is.null(mm)) 
        mm <- marksubset(mm, retain)
      result <- ppp(xoff, yoff, window = win, check = FALSE, 
                    marks = mm)
    }
    if (is.im(mu)) {
      P <- eval.im(mu/mumax)
      result <- rthin(result, P)
    }
    if (saveparents) {
      attr(result, "parents") <- parents
      attr(result, "parentid") <- parentid
      attr(result, "expand") <- expand
    }
    resultlist[[i]] <- result
  }
  result <- simulationresult(resultlist, nsim, drop)
  return(result)
}
getAvgSim<-function(sim){
  n<-length(sim)
  avgEnd<-0
  avgBr<-0
  avgDot<-0
  for(i in 1:n){
    avgEnd<-avgEnd+sim[[i]]$n/n
    avgBr<-avgBr+attr(sim[[i]],"parents")$n/n
    avgDot<-avgDot+(sim[[i]]$n)/(attr(sim[[i]],"parents")$n)/n
  }
  
  return(c(avgEnd,avgBr,mean(avgDot)))
}

