library(tidyverse)
library(ggplot2)
library(cowplot)
library(latex2exp)

# initialise results table
setwd("~/Desktop/Aalto/projpred-tsa/R")
csv_fname = "experiments/arma(2,1).csv"
df = read.table(file = csv_fname, sep = ",", header = TRUE)
fill_colour = "#999999"
true_colour = "#56B4E9"
ref_colour = "#E69F00"

p <- ggplot(df, aes(x=p_perp)) + 
  geom_histogram(fill=fill_colour) + 
  geom_vline(
    aes(xintercept=mean(p_true)), color=true_colour, linetype="dashed", size=1
  ) +
  geom_vline(
    aes(xintercept=5), color=ref_colour, linetype="dashed", size=1
  ) +
  ylab("Count") +
  xlab(TeX("$p^{perp}")) +
  annotate(x=mean(df$p_true), y=+Inf,label=TeX("$p^{true}$"), vjust=2, geom="label") +
  annotate(x=5, y=+Inf,label=TeX("$p^{*}$"), vjust=2, geom="label") +
  xlim(-0.5, 5.5) + 
  theme_bw()
q <- ggplot(df, aes(x=q_perp)) + 
  geom_histogram(fill=fill_colour) + 
  geom_vline(
    aes(xintercept=mean(q_true)), color=true_colour, linetype="dashed", size=1
  ) +
  geom_vline(
    aes(xintercept=5), color=ref_colour, linetype="dashed", size=1
  ) +
  ylab("Count") +
  xlab(TeX("$q^{perp}$")) +
  annotate(x=mean(df$q_true), y=+Inf,label=TeX("$q^{true}$"), vjust=2, geom="label") +
  annotate(x=5, y=+Inf,label=TeX("$q^{*}$"), vjust=2, geom="label") +
  xlim(-0.5, 5.5) +
  theme_bw()
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
# cowplot::plot_grid(p, q, P, Q, ncol = 2)
