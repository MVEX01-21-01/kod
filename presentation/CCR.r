library(ggplot2)
library(spatstat)
theme_set(theme_minimal())
source('plotting.r')

pC <- rMatClust(12, .07, 12)
pS <- rpoispp(12*12)
pR <- rMaternII(12*12, .1)
ppps <- list(pC, pS, pR)

pppplot2 <- function(ppp, title) {
  pppplot(ppp, parents=F) + ggtitle(title) +
    theme(plot.title=element_text(face = 'bold', size=16, hjust=0.5))
}
statplot <- function(ppps, stat, tf=function(r,y) y, ytxt=NULL, ...) {
  fvs <- lapply(ppps, stat, ...)
  xname <- attributes(fvs[[1]])$argu
  yname <- attributes(fvs[[1]])$valu
  if (is.null(ytxt)) ytxt <- attributes(fvs[[1]])$ylab
  ymax <- max(sapply(fvs, function(fv) max(tf(fv[, xname, drop=T], fv[, yname, drop=T]), tf(fv[, xname, drop=T], fv[, 'theo', drop=T]))))
  ymin <- min(sapply(fvs, function(fv) min(tf(fv[, xname, drop=T], fv[, yname, drop=T]), tf(fv[, xname, drop=T], fv[, 'theo', drop=T]))))
  
  lapply(fvs, function(fv) {
    df <- data.frame(r=fv[, xname, drop=T], y=fv[, yname, drop=T], theo=fv[, 'theo', drop=T])
    ggplot(df, aes(x=r, y=tf(r,y))) +
      geom_line(aes(y=tf(r,theo)), linetype=2) + geom_line(color='red') +
      ylab(ytxt) + ylim(ymin, ymax)
  })
}

grid.arrange(grobs=Map(pppplot2, ppps, title=c('Klustrad', 'CSR', 'Reguljär')), ncol=3)

grid.arrange(grobs=c(
  Map(pppplot2, ppps, title=c('Klustrad', 'CSR', 'Reguljär')),
  statplot(ppps, stat=Kest, correction='iso'),
  statplot(ppps, stat=Lest, correction='iso', tf=function(r,y) y - r, ytxt='L(r) - r')
  ), ncol=3)
