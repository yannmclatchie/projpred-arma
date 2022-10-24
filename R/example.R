# load required packages
# devtools::install_github("yannmclatchie/projpred@time_series")
library("forecast")
source("projpred_arma.R")

# simulate ARMA(3, 2) model data
series <- arima.sim(
  model = list(ar=c(0.9, -0.45, 0.25), ma=c(0.6, 0.45)), 
  n = 100
)
plot(series)

# use both procedures to compute submodel orders
projpred_lags <- projpred_arma(
  series, P_max = 5, Q_max = 5, stat = "elpd", silent = FALSE
)
print(projpred_lags)

# use auto.arima to extract submodel orders
auto_lags <- auto.arima(
  series, allowdrift=FALSE, max.d = 0, max.D = 0, max.P = 0, max.Q = 0
)
print(auto_lags)
