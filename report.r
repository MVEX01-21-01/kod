### This file should render the results of the report exactly,
### for reproducibility.

# Initialization ==========

if (!dir.exists('DATA_ENFS')) {
  stop('Can\'t find pattern subdir DATA_ENFS. Check your working directory.')
}

library(ggplot2)
theme_set(theme_minimal())
library(spatstat)

source('util.r')
source('repcluster.r')
source('multiGET.r')
data <- loaddata()

# Clean output directory
unlink('report_out', recursive=T)
dir.create('report_out')

# 3.0? - CSR test on patterns ==========
# TODO this is very ugly and unreadable, fix it
g <- groupplot(grouped(csrenvs, data, nsim=999))
ggsave('report_out/0.pdf', plot=g, width=6, height=6)

# 3.1 - Individual parameter fits ==========
fit.each.thomas   <- params.each(data$ppp, 'Thomas')
fit.each.matclust <- params.each(data$ppp, 'MatClust')
longparams <- function(pars, model, group) {
  # reshape the data into easily visualized format
  df <- data.frame(t(pars), model, group)
  varying <- colnames(df)[1:(ncol(df)-2)]
  df <- reshape(df, direction='long', varying=varying, v.names='value',
                timevar='param', times=varying)
  rownames(df) <- NULL
  df
}

data.params <- rbind(longparams(fit.each.thomas, 'Thomas', data$g),
                     longparams(fit.each.matclust, 'MatClust', data$g))
g <- ggplot(data.params, aes(group, value)) +
  geom_boxplot() + facet_wrap(~ model + param, scales='free')
ggsave('report_out/1.pdf', plot=g, width=5, height=5)

# 3.2 - Individual envelopes ==========
# TODO (just load them from disk)

# 3.3 - Group parameter fits ==========
# TODO what should we even do here? show the Kbar?

# 3.4 - Group envelopes ==========
# TODO (just load them from disk)
