# clean up environment
rm(list = ls())
gc(reset = TRUE)
setwd("/Users/yannmclatchie/Desktop/Aalto/projpred-arma/R/experiments")

library(tidyverse)
library(brms)
library(simts)
library(forecast)

# define experimental parameters
num_sims <- 100
sim_length <- 500
P_max <- 5
Q_max <- P_max

# define model
# model <- list(ar = c(0.5))
model <- list(ma = c(0.5, 0.3))
# model <- list(ar=c(-0.2,-0.3,0.5), ma=c(0.4,.3))
# model <- list(ar=c(0.2, -0.1), ma=c(0.6))
# model <-list(ar = c(0.9), ma = c(2, 0.5))
P_true <- length(model$ar)
Q_true <- length(model$ma)

# initialise results table
csv_fname = paste0("auto-arma(",P_true,",",Q_true,").csv")
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c('P_star', 'Q_star', 'P_true', 'Q_true', 'P_perp', 'Q_perp')
##
## WARNING: this will rewrite any existing table
##
write.table(df, file = csv_fname, sep = ",") 

# run tests
num_steps <- num_sims
ii <- 0
for (sim in 1:num_sims){
  cat(paste0(round(ii / num_steps * 100), '% completed'))
  # simulated data set
  data <- arima.sim(model = model, n = sim_length, n.start = 100)
  # get auto.arima predictions
  auto.arima.fit <- auto.arima(
    data,
    max.p = P_max,
    max.q = Q_max,
    max.P = 0,
    max.Q = 0,
    max.order = 10,
    max.d = 0,
    max.D = 0
  )
  P_perp <- length(auto.arima.fit$model$phi)
  Q_perp <- length(auto.arima.fit$model$theta)
  # append results to experiments table
  row <- data.frame(
    t(
      c(
        P_max,
        Q_max,
        P_true, 
        Q_true, 
        P_perp, 
        Q_perp
      )
    )
  )
  write.table(row, file = csv_fname, sep = ",",
              append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  if (ii == num_steps) cat(': Done')
  else cat('\014')
  ii <- ii + 1
}

library(ggplot2)
library(cowplot)
library(latex2exp)

df = read.table(file = csv_fname, sep = ",", header = TRUE)
fill_colour = "#999999"
true_colour = "#56B4E9"
ref_colour = "#E69F00"

P <- ggplot(df, aes(x=P_perp)) + 
  geom_histogram(fill=fill_colour) + 
  geom_vline(
    aes(xintercept=mean(P_true)), color=true_colour, linetype="dashed", size=1
  ) +
  geom_vline(
    aes(xintercept=mean(P_star)), color=ref_colour, linetype="dashed", size=1
  ) +
  ylab("Count") +
  #xlab(TeX("$P^{perp}")) +
  xlab(TeX("$p^{perp}")) +
  annotate(x=mean(df$P_true), y=+Inf,label=TeX("$p^{true}$"), vjust=2, geom="label") +
  annotate(x=5, y=+Inf,label=TeX("$p^{*}$"), vjust=2, geom="label") +
  #annotate(x=mean(df$P_true), y=+Inf,label=TeX("$P^{true}$"), vjust=2, geom="label") +
  #annotate(x=mean(df$P_star), y=+Inf,label=TeX("$P^{*}$"), vjust=2, geom="label") +
  xlim(-0.5, mean(df$P_star)+0.5) + 
  theme_bw()
Q <- ggplot(df, aes(x=Q_perp)) + 
  geom_histogram(fill=fill_colour) + 
  geom_vline(
    aes(xintercept=mean(Q_true)), color=true_colour, linetype="dashed", size=1
  ) +
  geom_vline(
    aes(xintercept=mean(Q_star)), color=ref_colour, linetype="dashed", size=1
  ) +
  ylab("Count") +
  #xlab(TeX("$Q^{perp}$")) +
  xlab(TeX("$q^{perp}")) +
  annotate(x=mean(df$Q_true), y=+Inf,label=TeX("$q^{true}$"), vjust=2, geom="label") +
  annotate(x=5, y=+Inf,label=TeX("$q^{*}$"), vjust=2, geom="label") +
  #annotate(x=mean(df$Q_true), y=+Inf,label=TeX("$Q^{true}$"), vjust=2, geom="label") +
  #annotate(x=mean(df$Q_star), y=+Inf,label=TeX("$Q^{*}$"), vjust=2, geom="label") +
  xlim(-0.5, mean(df$Q_star)+0.5) +
  theme_bw()
cowplot::plot_grid(P, Q, ncol = 1)
