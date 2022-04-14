# clean up environment
rm(list = ls())
gc(reset = TRUE)
# build and install the local version of projpred
setwd("/Users/yannmclatchie/Desktop/Aalto/projpred")
getwd()
devtools::build()
devtools::install()

library(projpred)
library(tidyverse)
library(brms)
library(simts)
library(latex2exp)
setwd("~/Desktop/Aalto/projpred-tsa/R")
source("arma_projpred.R")

# define experimental parameters
sim_length <- 300
P <- 6

# simulate and plot data from model 
model <- list(ar=c(-0.8, 0.5, 0.3, -0.1, -0.04, 0.01))
model <- list(ar=c(0.9, -0.45, 0.25, -0.08, -0.04, 0.01))
data <- arima.sim(model = model, n = sim_length, n.start = 100)
plot(data, type="l")

# plot "true" model parameter values
df <- data.frame(model$ar)
df$lag <- 1:length(df$model.ar)
df %>%
  ggplot(aes(x=lag, y=model.ar)) +
  geom_point(size = 1, colour = "black") + 
  geom_hline(yintercept=0)+ 
  geom_segment( aes(x=lag, xend=lag, y=0, yend=model.ar)) +
  labs(x=TeX("Lag, $p$"), y=TeX("Value of parameter, $\\theta_p$")) +
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
# build design matrix
X <- as.data.frame(embed(data, P+1)) %>%
  rev() %>%
  `colnames<-`(sprintf("x.t%s", 0:P))
# build AR formula
if (P == 0) {
  ar_formula <- formula("x.t0 ~ 1")
} else {
  lagged_covs <- colnames(X)[colnames(X) != "x.t0"]
  ar_formula <- formula(
    paste("x.t0 ~ 1 + ", paste0(lagged_covs, collapse = " + "))
  )
}


# fit sparse reference model
ar_sparse <- brms::brm(
  ar_formula,
  data = X,
  family = gaussian,
  prior = brms::prior(horseshoe()),
  silent = 2,
  backend = "rstan",
  open_progress = FALSE,
  refresh = 0
)
ref_loo <- max(sparse_vs$summary$elpd)
# apply projpred
sparse_ref <- projpred::get_refmodel(ar_sparse)
sparse_vs <- projpred::varsel(sparse_ref, method = "arma", verbose = FALSE)
plot(sparse_vs)

p1 <- df %>%
  ggplot(aes(x=lag, y=model.ar)) +
  geom_point(size = 2, colour = "black") + 
  geom_hline(yintercept=0)+ 
  geom_segment( aes(x=lag, xend=lag, y=0, yend=model.ar)) +
  labs(x=TeX("Lag, $p$"), y=TeX("Value of parameter, $\\theta_p$")) +
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.ticks.x = element_blank())+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  expand_limits(y=c(-0.4, 0.8))+
  expand_limits(x=c(0, 6))+
  scale_y_continuous(breaks=seq(-0.4, 0.8, 0.2))
p2 <- ggplot(data=sparse_vs$summary, aes(x=size, y=elpd)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=ref_loo, color="red", linetype = "dashed")+
  geom_errorbar(aes(ymin=elpd-se, ymax=elpd+se), width=0)+
  labs(x=TeX("Model size"), y=TeX("ELPD")) +
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.y = element_line(color="black", size = 0.5))+
  theme(axis.line.x = element_line(color="black", size = 0.5))
cowplot::plot_grid(p1, p2, ncol = 1)

# fit wide reference model
ar_wide <- brms::brm(
  ar_formula,
  data = X,
  family = gaussian,
  prior = brms::prior_string("normal(0,10)"),
  silent = 2,
  backend = "rstan",
  open_progress = FALSE,
  refresh = 0
)
# apply projpred
wide_ref <- projpred::get_refmodel(ar_wide)
wide_vs <- projpred::varsel(wide_ref, method = "arma", verbose = FALSE)
plot(wide_vs)
