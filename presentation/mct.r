obsppp <- data$ppp[[3]]
fakeppp <- data$ppp[[5]]
testM <- testM <- kppm(obsppp, cluster='Thomas')
sims <- c(simulate(testM, nsim=19))
grid.arrange(grobs=Map(pppplot, c(list(fakeppp), sims[1:8]), parents=F))

csetfv <- function(fvs, obs) {
  xname <- attributes(obs)$argu
  yname <- attributes(obs)$valu
  create_curve_set(list(
    r=fvs[[1]][, xname, drop=T],
    obs=obs[, yname, drop=T],
    sim_m=cbind(sapply(fvs, function(fv) fv[, yname, drop=T]))))
}
curveG <- function(fv, txt) {
  xname <- attributes(fv)$argu
  yname <- attributes(fv)$valu
  df <- data.frame(r=fv[, xname, drop=T], G=fv[, yname, drop=T])
  ggplot(df, aes(x=r, y=G)) +
    geom_line() + geom_text(data=data.frame(txt=txt),x=55,y=0.25,label=txt,size=8)
}

Gs <- lapply(sims, Gest, correction='km')
cset <- csetfv(Gs, Gest(fakeppp, correction='km'))
plot(cset) + theme(legend.position='none') + 

df <- as.data.frame(cbind(x=cset$r, cset$funcs)) %>%
  gather(key='id', value='y', -x)
erl <- forder(cset)
ggplot(filter(df, id!='obs'), aes(x=x, y=y, group=id, color=rep(erl[-1], each=513))) +
  geom_line() +
  geom_line(data=filter(df, id=='obs'), aes(group='obs'), color='white', size=2.5) +
  geom_line(data=filter(df, id=='obs'), aes(group='obs', color=erl[1]), size=1) +
  labs(x='r', y='G(r)', color='ERL')

plot(global_envelope_test(cset,alpha=1/(length(sims)+1)))


####

Klimits <- list(
  xmin = max(sapply(K, function(fv) min(fv$r))),
  xmax = min(sapply(K, function(fv) max(fv$r))),
  ymin = min(sapply(K, function(fv) min(fv$loiso, fv$lotheo, na.rm=T))),
  ymax = max(sapply(K, function(fv) max(fv$hiiso, fv$hitheo, na.rm=T)))
)
Kgs <- Map(function(fv, th, mc) {
  df <- as.data.frame(fv)
  ggplot(df, aes(r, pooliso)) +
    geom_line(aes(colour='empirical')) +
    geom_ribbon(aes(ymin=loiso,ymax=hiiso),alpha=0.1, colour=c("#e0e0e0")) +
    geom_line(aes(,pooltheo), linetype = 'dashed', colour=c("#c0c0c0")) +
    geom_line(aes(,th$fit$fit, colour='Thomas'), linetype='dashed') +
    geom_line(aes(,mc$fit$fit, colour='MatClust'), linetype='dashed') +
    coord_cartesian(xlim=c(Klimits$xmin, Klimits$xmax), ylim=c(Klimits$ymin, Klimits$ymax)) +
    scale_colour_manual(values=c(c("#444444"),c("#F7951F"),c("#0154A6"))) +
    xlab(expression(italic(r))) + ylab(expression(italic(K(r)))) +
    theme(legend.position='bottom', legend.title=element_blank())
}, K, fit.thomas.K, fit.matclust.K)
Kaux <- auxgrobs(Kgs[[1]])
g <- combi.grouped(lapply(Kgs, function(g) list(g + labs(x=NULL, y=NULL) + theme(legend.position='none'))), aux=Kaux)
