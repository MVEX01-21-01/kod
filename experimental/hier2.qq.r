# 3.7 - hier2.R ==========
# Alexander's qqplot
H2 <- new.env()
source('experimental/hier2.R', local=H2)
aqq.genq <- function(data, cluster, group, estfn) {
  ys <- do.call(estfn, data)
  data.frame(
    cluster=cluster,
    group=group,
    value=ys
  )
}
dfs.aqq <- Map(aqq.genq,
               data=rep(list(list(H2$data_moderate, H2$data_moderate_b), list(H2$data_normal, H2$data_normal_b)), length.out=4),
               cluster=rep(c('Thomas','MatClust'),each=2),
               group=rep(c('MODERATE','NORMAL'),length.out=4),
               estfn=rep(list(
                 H2$thomas.scale.est2,
                 function(...) Filter(function(x) x>0, H2$mat.scale.est2(...))
               ), each=2)
)
distrs.aqq <- {
  sdM <- H2$thomas.scale.est(H2$data_moderate, H2$data_moderate_b)
  sdN <- H2$thomas.scale.est(H2$data_normal, H2$data_normal_b)
  RM <- H2$mat.scale.est(H2$data_moderate, H2$data_moderate_b)
  RN <- H2$mat.scale.est(H2$data_normal, H2$data_normal_b)
  list(
    function(p) sdM*sqrt(-2*log(1-p)),
    function(p) sdN*sqrt(-2*log(1-p)),
    function(p) qunif(p, max=RM),
    function(p) qunif(p, max=RN)
  )
}
g <- grid.arrange(grobs=Map(function(data, distr) 
  ggplot(data, aes(sample=value)) +
    stat_qq(distribution=distr) + stat_qq_line(distribution=distr) + facet_wrap(~ cluster + group),
  dfs.aqq, distr=distrs.aqq))
