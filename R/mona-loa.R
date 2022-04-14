# clean up environment
rm(list = ls())
gc(reset = TRUE)
# build and install the local version of projpred
setwd("/Users/yannmclatchie/Desktop/Aalto/projpred")
getwd()
devtools::build()
devtools::install()

# remotes::install_github("asael697/bayesforecast", dependencies = TRUE)

library(projpred)
library(bayesforecast)
library(tidyverse)
library(brms)
library(simts)
setwd("~/Desktop/Aalto/projpred-tsa/R")
source("arma_projpred.R")

library(ggfortify)
library(cowplot)

## 1. ARMA model 

# load and plot data
data("LakeHuron")
# clean data
data <- ts(LakeHuron)

p2 <- autoplot(acf(data, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))

p3 <- autoplot(pacf(data, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  ylab("PACF") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))

p1 <- ggplot(data.frame(Y=as.matrix(data), date=time(data)), aes(y=Y, x=date)) +
  geom_line() + 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  xlab("Date") +
  theme(legend.position = 'bottom')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )

cowplot::plot_grid(p1, p2, p3, ncol = 1)

data_diff <- data %>%
  diff()
acf(data_diff)
pacf(data_diff)
plot(data_diff)


p1_diff <- ggplot(data.frame(Y=as.matrix(data_diff), date=time(data_diff)), aes(y=Y, x=date)) +
  geom_line() + 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  xlab("Date") +
  theme(legend.position = 'bottom')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )
p2_diff <- autoplot(acf(data_diff, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
p3_diff <- autoplot(pacf(data_diff, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  ylab("PACF") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
cowplot::plot_grid(p1, p1_diff, p2, p2_diff, p3, p3_diff, ncol = 2)



P_max <- 4
Q_max <- 4
restricted_lags <- proj_arma(
  data_diff, P_max=P_max, Q_max=Q_max, pct = 0.1, stat = "elpd", silent=TRUE
)
message(restricted_lags)

# fit reference model
reference_model <- brms::brm(
  x ~ arma(
    p = P_max, 
    q = Q_max
  ),
  data = data_diff
)
reference_model <- add_criterion(reference_model, c("loo", "waic"))
loo(reference_model)
# fit restricted model
restricted_model <- brms::brm(
  x ~ arma(
    p = 1, 
    q = 1
    ),
  data = data_diff
  )
restricted_model <- add_criterion(restricted_model, c("loo", "waic"))
loo(restricted_model)
# plot pareto k values
plot(loo(reference_model), diagnostic = "k", label_points = TRUE)
plot(loo(restricted_model), diagnostic = "k", label_points = TRUE)
# compare the two
loo_compare(reference_model, restricted_model)
waic(reference_model, restricted_model)
loo(reference_model, restricted_model)

## 2. SARMA model

data("co2")
plot(co2)
# clean data
data <- ts(co2) %>% 
  diff(lag=12) %>%
  diff()
acf(data)
pacf(data)
plot(data)

p1 <- ggplot(data.frame(Y=as.matrix(co2), date=time(co2)), aes(y=Y, x=date)) +
  geom_line() + 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  xlab("Date") +
  theme(legend.position = 'bottom')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )
p2 <- autoplot(acf(co2, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
p3 <- autoplot(pacf(co2, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  ylab("PACF") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
p1_diff <- ggplot(data.frame(Y=as.matrix(data), date=time(data)), aes(y=Y, x=date)) +
  geom_line() + 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  xlab("Date") +
  theme(legend.position = 'bottom')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )
p2_diff <- autoplot(acf(data, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
p3_diff <- autoplot(pacf(data, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') +
  geom_hline(yintercept=0)+ 
  scale_color_manual(values=c('blue','magenta','red','green'), 
                     breaks=c('Profitability', 'Growth', 'Safety','Payout'))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  ylab("PACF") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
cowplot::plot_grid(p1, p1_diff, p2, p2_diff, p3, p3_diff, ncol = 2)

p_max <- 4
q_max <- p_max
P_max <- 2
Q_max <- P_max

seasonality <- 12
seasonal_x <- data[seq(1, length(data), seasonality)]
non_seasonal_x <- data
# run projpred sarma procedure
seasonal_restricted_lags <- proj_arma(
  seasonal_x, P_max=P_max, Q_max=Q_max, pct = 0.1, stat = "elpd", silent=TRUE
)
non_seasonal_restricted_lags <- proj_arma(
  non_seasonal_x, P_max=p_max, Q_max=q_max, pct = 0.1, stat = "elpd", silent=TRUE
)
message(seasonal_restricted_lags, non_seasonal_restricted_lags)

reference = stan_sarima(
  ts = data, 
  order = c(p_max, 0, q_max), 
  seasonal = c(P_max, 0, Q_max)
)
plot(loo(reference), diagnostic = "k", label_points = TRUE)
print(waic(reference))


restricted = stan_sarima(
  ts = data, 
  order = c(non_seasonal_restricted_lags$P_perp, 0, non_seasonal_restricted_lags$Q_perp), 
  seasonal = c(seasonal_restricted_lags$P_perp, 0, seasonal_restricted_lags$Q_perp)
)
plot(loo(restricted), diagnostic = "k", label_points = TRUE)
print(waic(restricted))

loo_compare(loo(reference), loo(restricted))

loo(reference)
loo(restricted)


