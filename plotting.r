library(ggplot2)
library(gridExtra)
library(grid)

# General ====
# Relabeller
relabel <- function(main, xlab, ylab, legend) {
  if (!is.null(legend)) {
    lheight <- sum(legend$height)
    xheight <- sum(xlab$height)
    ywidth <- sum(ylab$width)
    
    grid.arrange(main, xlab, ylab, legend,
                 layout_matrix=rbind(c(3,1),c(NA,2),c(NA,4)),
                 heights=unit.c(unit(1, "npc") - lheight - xheight, xheight, lheight),
                 widths=unit.c(ywidth, unit(1, "npc") - ywidth))
  } else {
    xheight <- sum(xlab$height)
    ywidth <- sum(ylab$width)
    
    grid.arrange(main, xlab, ylab,
                 layout_matrix=rbind(c(3,1),c(NA,2)),
                 heights=unit.c(unit(1, "npc") - xheight, xheight),
                 widths=unit.c(ywidth, unit(1, "npc") - ywidth))
  }
}

# Extract labels and legend
auxgrobs <- function(g) {
  grs <- ggplotGrob(g)$grobs
  legendq <- which(sapply(grs, function(x) x$name) == "guide-box")
  list(
    legend = if (length(legendq) > 0) grs[[legendq]] else NULL,
    xlab = grs[[which(sapply(grs, function(x) grepl('title.x', x$name, fixed=T)))]],
    ylab = grs[[which(sapply(grs, function(x) grepl('title.y', x$name, fixed=T)))]]
  )
}

# Combine
combi <- function(gs, aux=NULL, ...) {
  if (is.null(aux)) { aux <- auxgrobs(gs[[1]]) }
  g <- arrangeGrob(grobs=gs, ...)
  relabel(g, aux$xlab, aux$ylab, aux$legend)
}

# Combine, grouped
combi.grouped <- function(groups, aux=NULL, ...) {
  if (is.null(aux)) { aux <- auxgrobs(groups[[1]]) }
  g <- arrangeGrob(grobs=Map(function(group, name) {
    arrangeGrob(grobs=group, top=name, ...)
  }, groups, names(groups)), ncol=length(groups))
  relabel(g, aux$xlab, aux$ylab, aux$legend)
}

# Envelopes ====
# Extract a raw envelope grob
envelopegrob <- function(envelope, tag=NULL, noaxes=F) {
  g <- plot(envelope) +
    switch(noaxes+1, theme(legend.position='none'), theme(legend.position='none', axis.text.x=element_blank(), axis.text.y=element_blank())) +
    labs(title=NULL, x=NULL, y=NULL)
  if (!is.null(tag)) {
    g <- g +
      labs(tag=paste0('#',tag)) +
      theme(plot.tag.position = c(0.05, 0.95), plot.tag=theme_get()$axis.text)
  }
  g
}

# Extract lims from envelope
envxlim <- function(envelope) {
  c(envelope$r[1], envelope$r[length(envelope$r)])
}
envylim <- function(envelope) {
  c(min(envelope$obs, envelope$lo), max(envelope$obs, envelope$hi))
}

# Combine envelopes
plot.envs.single <- function(envelopes, noaxes=T, ...) {
  # TODO common limits?
  aux <- auxgrobs(plot(envelopes[[1]]))
  combi(Map(envelopegrob, envelopes, tag=names(envelopes), noaxes=noaxes), aux, ...)
}

# Combine envelopes grouped
plot.envs.grouped <- function(groups, noaxes=T, ...) {
  aux <- auxgrobs(plot(groups[[1]][[1]]))
  combi.grouped(lapply(groups, function(group) Map(envelopegrob, group, tag=names(group), noaxes=noaxes)), aux, ...)
}

# Point patterns ====
pppplot <- function(X, tag=NULL, parents=T, tree=F) {
  Xm <- NULL
  Pm <- NULL
  if (!is.null(attr(X, 'parents')) && parents) {
    P <- attr(X, 'parents')
    if (tree) {
      Xm <- factor(attr(X, 'parentid'))
      Pm <- factor(1:npoints(P))
    }
  } else {
    P <- list(x=numeric(0),y=numeric(0))
  }
  if (!is.null(tag)) tag <- paste0('#',tag)

  xlim <- X$window$xrange
  ylim <- X$window$yrange
  ggplot() + geom_point(aes(x=X$x, y=X$y, color=Xm), shape=1, stroke=0.35) +
    geom_point(aes(x=P$x, y=P$y, color=Pm), size=3, shape=4, stroke=0.35) +
    labs(title=tag, tag=NULL, x=NULL, y=NULL) + coord_fixed(xlim=xlim, ylim=ylim, expand=F, clip='off') +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(color='gray', fill=NA, size=0.5), plot.title=theme_get()$axis.text, legend.position='none')
}
