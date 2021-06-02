# ERL example

L <- 9
N <- 6
maxv <- 10

X <- 1:L
Y <- t(replicate(L, sample.int(maxv, size=N)))
cs <- create_curve_set(list(r=X, obs=Y))

cp <- data.frame(x=rep.int(X, N), y=c(Y), series=factor(rep(1:N, each=L)))
g <- ggplot(cp,aes(x=x,y=y)) +
  geom_line(aes(color=series)) + geom_point(aes(color=series), shape=15, size=2) +
  scale_x_continuous(breaks=X) + scale_y_continuous(breaks=1:maxv) +
  labs(x='j', y='S_ij', color='i') +
  theme(panel.grid.major = element_line(colour='grey70'))

ggsave('Aerl_sampledata.pdf', plot=g, width=5, height=2.1)

pwer <- function(v) {
  ranks <- rank(v, ties.method = "average")
  pmin(ranks, length(v) + 1 - ranks)
}
cpw <- apply(Y, 1, pwer)
ckw <- t(apply(cpw, 1, sort))

forder(cs, measure='erl')
