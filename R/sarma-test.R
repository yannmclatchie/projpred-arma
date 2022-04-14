# clean up environment
rm(list = ls())
gc(reset = TRUE)
# build and install the local version of projpred
setwd("/Users/yannmclatchie/Desktop/Aalto/projpred")
getwd()
devtools::build()
devtools::install()

devtools::install_github("SMAC-Group/simts")
library(simts)
library(projpred)
library(tidyverse)
library(brms)
setwd("~/Desktop/Aalto/projpred-tsa/R")
source("arma_projpred.R")

# define model
non_seasonal_model <- list(ar = c(0.5))
seasonal_model <- list(ar = c(0.9), ma = c(2, 0.5))
seasonality <- 12
p_true <- length(non_seasonal_model$ar)
q_true <- length(non_seasonal_model$ma)
P_true <- length(seasonal_model$ar)
Q_true <- length(seasonal_model$ma)

# define model
model <- SARMA(
  ar=non_seasonal_model$ar, 
  ma=non_seasonal_model$ma, 
  sar=seasonal_model$ar, 
  sma=seasonal_model$ma, 
  s=seasonality
)

# define experimental parameters
num_sims <- 100
N <- 500
sim_length <- N * seasonality
p_max <- 5
q_max <- q_max
P_max <- 3
Q_max <- P_max

# initialise results table
csv_fname = paste0(
  "experiments/arma(",
  p_true,
  ",",
  q_true,
  ")x("
  ,P_true,
  ",",
  Q_true,
  ")_",
  seasonality,
  ".csv"
)
df <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(df) <- c(
  'p_star', 
  'q_star', 
  'p_true', 
  'q_true', 
  'p_perp', 
  'q_perp', 
  'P_star', 
  'Q_star', 
  'P_true', 
  'Q_true', 
  'P_perp', 
  'Q_perp'
)
##
## WARNING: this will rewrite any existing table
##
# write.table(df, file = csv_fname, sep = ",") 

# run experiment
num_steps <- num_sims
ii <- 0
for (sim in 1:num_sims){
    cat(paste0(round(ii / num_steps * 100), '% completed'))
    # simulated data set
    x = gen_gts(sim_length, model)[,1]
    seasonal_x <- x[seq(1, length(x), seasonality)]
    non_seasonal_x <- x[1:500]
    # run projpred arma procedure
    seasonal_restricted_lags <- proj_arma(
      seasonal_x, P_max=P_max, Q_max=Q_max, pct = 0.1, stat = "elpd", silent=TRUE
    )
    non_seasonal_restricted_lags <- proj_arma(
      non_seasonal_x, P_max=p_max, Q_max=q_max, pct = 0.1, stat = "elpd", silent=TRUE
    )
    # append results to experiments table
    row <- data.frame(
      t(
        c(
          P_max, 
          Q_max, 
          p_true, 
          q_true, 
          non_seasonal_restricted_lags$P_perp, 
          non_seasonal_restricted_lags$Q_perp,
          P_max, 
          Q_max, 
          P_true, 
          Q_true, 
          seasonal_restricted_lags$P_perp, 
          seasonal_restricted_lags$Q_perp
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
