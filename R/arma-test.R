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
setwd("~/Desktop/Aalto/projpred-tsa/R")
source("arma_projpred.R")

# define experimental parameters
num_sims <- 100
sim_length <- 500
P_max <- 5
Q_max <- P_max

# define model
# model <- list(ar = c(0.5))
# model <- list(ma = c(0.5, 0.3))
# model <- list(ar=c(-0.2,-0.3,0.5), ma=c(0.4,.3))
model <- list(ar=c(-0.5, 0.9), ma=c(0.6))
# model <-list(ar = c(0.9), ma = c(2, 0.5))
P_true <- length(model$ar)
Q_true <- length(model$ma)

# initialise results table
csv_fname = paste0("experiments/arma(",P_true,",",Q_true,").csv")
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c('P_star', 'Q_star', 'P_true', 'Q_true', 'P_perp', 'Q_perp')
##
## WARNING: this will rewrite any existing table
##
# write.table(df, file = csv_fname, sep = ",") 

# run tests
num_steps <- num_sims
ii <- 0
for (sim in 1:num_sims){
    cat(paste0(round(ii / num_steps * 100), '% completed'))
    # simulated data set
    data <- arima.sim(model = model, n = sim_length, n.start = 100)
    restricted_lags <- proj_arma(
      data, P_max=P_max, Q_max=Q_max, pct = 0.1, stat = "elpd", silent=TRUE
    )
    # get auto.arima predictions
    auto.arima.fit <- auto.arima(data)
    # append results to experiments table
    row <- data.frame(
      t(
        c(
          P_max, 
          Q_max, 
          P_true, 
          Q_true, 
          restricted_lags$P_perp, 
          restricted_lags$Q_perp
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
