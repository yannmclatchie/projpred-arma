# Kullback-Leibler Projections for Automatic Order Identification of Seasonal ARMA Models

This repository implements Algorithm 1 of the paper _Bayesian order identification of ARMA models with projection predictive inference_ by McLatchie et al. ([2022](https://arxiv.org/abs/2208.14824)).

## Algorithm

**Inputs:** Time series data $\{y_t\}$, reference model $P^\ast$, reference model $Q^\ast$, [projpred model selection heuristic](https://mc-stan.org/projpred/reference/suggest_size.html) parameters `alpha` and `pct`

1. Perform Ljung-Box test for stationarity to data $\{y_t\}$ and print warning message if failed
2. Fit a linear reference model (AR) to these data observations $\{y_t\}$ with `BRMS` and lag parameter $P^\ast$
3. Apply the model selection heuristic to this reference AR model to get some restricted AR lag value $P^\perp$
4. Fit a linear restricted model (AR) to $\{y_t\}$ with `BRMS` and this new restricted lag parameter $P^\perp$
5. Extract the residuals from this AR restricted model, denote these residuals $\{\epsilon_t\}$
6. Fit a linear reference model (MA) to these residuals $\{\epsilon_t\}$ with `BRMS` and lag parameter $Q^\ast$
7. Apply the model selection heuristic to this reference MA model to get some restricted MA lag value $Q^\perp$

**Returns:** The restricted lag values $P^\perp$ and $Q^\perp$

## Development

This algorithm requires a custom search heuristic, detailed in Figure 1 by McLatchie et al. ([2022](https://arxiv.org/abs/2208.14824)), which is available in a personal fork of `projpred` available for download through

```R
devtools::install_github("yannmclatchie/projpred@time_series)
```