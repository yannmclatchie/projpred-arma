require(devtools)
require(brms)

# load local version of projpred
setwd("~/Desktop/Aalto/projpred-arma")
if("projpred" %in% (.packages())){
  detach("package:projpred", unload=TRUE) 
}
load_all("./projpred/")

# define default arma prior
arma_prior <- c(
  prior_string("normal(0, 0.5)", class = "b"),
  # prior_string("student_t(7, 7, 3)", class = "sigma"),
  prior_string("student_t(6, 10, 3)", class = "sigma"),
  prior_string("student_t(6, 0, 2.5)", class = "Intercept")
)

.fit_ar_model <- function(data, P){
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
  # fit AR model in BRMS
  ar_fit <- brm(
    formula = ar_formula,
    data = X,
    family = gaussian,
    prior = arma_prior,
    silent = 2,
    backend = "rstan",
    open_progress = FALSE,
    refresh = 0
  )
  return(ar_fit)
}

.fit_ma_model <- function(data, Q){
  # build design matrix of residuals
  E <- as.data.frame(embed(data, Q+1)) %>%
    rev() %>%
    `colnames<-`(sprintf("eps.t%s", 0:Q))
  # build MA formula
  if (Q == 0) {
    ma_formula <- formula("eps.t0 ~ 1")
  } else {
    lagged_covs <- colnames(E)[colnames(E) != "eps.t0"]
    ma_formula <-
      formula(paste("eps.t0 ~ 1 + ", paste0(lagged_covs, collapse = " + ")))
  }
  # fit MA model to residuals
  ma_fit <- brm(
    ma_formula,
    data = E,
    family = gaussian,
    prior = arma_prior,
    silent = 2,
    backend = "rstan",
    open_progress = FALSE,
    refresh = 0
  )
  return(ma_fit)
}

.get_ar_refmodel <- function(data, P_max = 6){
  # fit full model in BRMS
  brms_ar <- .fit_ar_model(data, P_max)
  # return refmodel
  ar_ref <- get_refmodel(brms_ar)
  return(ar_ref)
}

.get_ma_refmodel <- function(eps, Q_max = 6){
  # fit MA model to residuals
  brma_ma <- .fit_ma_model(eps, Q_max)
  # return refmodel
  ma_ref <- get_refmodel(brma_ma)
  return(ma_ref)
}

.get_restricted_lag <- function(
  refmodel, data, alpha = 0.32, pct = 0, cv = FALSE, stat = "elpd"
){
  # perform latent projection predictive variable selection on model
  if (!cv) {
    vs <- varsel(refmodel, method = "arma", verbose = FALSE)
  }
  else {
    vs <- cv_varsel(refmodel, method = "arma", verbose = FALSE)
  }
  # use heuristic to chose informative parameters in submodel
  suggested_size <- suggest_size(vs, alpha = alpha, pct = pct, stat = stat)
  # extract the largest lag from the restricted covariate space
  lag_perp <- length(
    unlist(as.list(solution_terms(vs))[0:suggested_size])
  )
  return(lag_perp)
}

# primary function
proj_arma <- function(
  data,
  P_max = 6,
  Q_max = 6,
  p_crit = 0.01,
  alpha = 0.32,
  pct = 0,
  cv = FALSE,
  stat = "elpd",
  silent = FALSE
){
  start.time <- Sys.time()
  
  # Ljung-Box test for stationarity
  p_val <- Box.test(data, type="Ljung-Box")$p.value
  if ((p_val < p_crit) && (!silent)) {
    message(paste0("Ljung-Box test p-value = ", round(p_val, 4)))
    message("Your time series may not be stationary, proceed with caution.\n")
  }
  
  # apply projpred to AR component
  if(!silent){message("Applying projpred to the AR component ... ")}
  ar_refmodel <- .get_ar_refmodel(data, P_max = P_max)
  # perform variable selection on AR component
  if (!cv) {
    vs <- varsel(ar_refmodel, method = "arma", verbose = FALSE)
  }
  else {
    vs <- cv_varsel(ar_refmodel, method = "arma", verbose = FALSE)
  }
  # use heuristic to chose informative parameters in submodel
  suggested_size <- suggest_size(vs, alpha = alpha, pct = pct, stat = stat)
  # extract the largest lag from the restricted covariate space
  P_perp <- length(
    unlist(as.list(solution_terms(vs))[0:suggested_size])
  )
  if(!silent){message(paste0("Done!\nP_perp = ", P_perp))}
  
  # extract residuals from projected AR component
  projected_ar <- project(vs, nterms = suggested_size)
  proj_ar_pred <- colMeans(proj_linpred(projected_ar)$pred)
  res <- proj_ar_pred - head(data, length(proj_ar_pred))
  
  # apply projpred to MA component
  if(!silent){message("Applying projpred to the MA component ... ")}
  ma_refmodel <- .get_ma_refmodel(res, Q_max = Q_max)
  Q_perp <- .get_restricted_lag(
    ma_refmodel, data, alpha = alpha, pct = pct, cv = cv, stat = stat
  )
  if(!silent){message(paste0("Done!\nQ_perp = ", Q_perp))}
  if(!silent){
    message(
      paste0(
        "Time elapsed: ~",
        as.integer(as.numeric(Sys.time() - start.time, units = "mins")),
        " minutes"
      )
    )
  }
  
  # return restricted lags
  return(list(P_perp = P_perp, Q_perp = Q_perp))
}
