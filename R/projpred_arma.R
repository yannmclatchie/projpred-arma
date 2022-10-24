library("projpred")
library("brms")
library("dplyr")

arma_prior <- c(
  prior_string("normal(0, 0.5)", class = "b"),
  prior_string("student_t(7, 0, 1)", class = "sigma"),
  prior_string("student_t(6, 0, 2.5)", class = "Intercept")
)

.fit_ar_model <- function(data, P){
  # build design matrix
  X <- as.data.frame(embed(data, P+1)) %>%
    rev() %>%
    `colnames<-`(sprintf("x.t%s", 0:P))
  # build AR formula
  if (P == 0) {
    ar_formula <- "x.t0 ~ 1"
  } else {
    lagged_covs <- colnames(X)[colnames(X) != "x.t0"]
    ar_formula <- paste(
      "x.t0 ~ 1 + ", paste0(lagged_covs, collapse = " + ")
    )
  }
  # fit AR model in BRMS
  ar_fit <- brms::brm(
    formula = ar_formula,
    data = X,
    prior = arma_prior,
    silent = 2,
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
  ma_fit <- brms::brm(
    ma_formula,
    data = E,
    prior = arma_prior,
    silent = 2,
    open_progress = FALSE,
    refresh = 0
  )
  return(ma_fit)
}

.perform_search <- function(refmodel, cv = FALSE){
  # perform latent projection predictive variable selection on model
  if (!cv) {
    vs <- projpred::varsel(refmodel, method = "arma", verbose = FALSE)
  }
  else {
    vs <- projpred::cv_varsel(refmodel, method = "arma", verbose = FALSE)
  }
  return(vs)
}

projpred_arma <- function(
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
  ar_refmodel <- .fit_ar_model(data, P = P_max)
  ar_search <- .perform_search(ar_refmodel)
  # use heuristic to chose informative parameters in submodel
  suggested_size <- projpred::suggest_size(
    ar_search, alpha = alpha, pct = pct, stat = stat
  )
  # extract the largest lag from the restricted covariate space
  P_perp <- length(
    unlist(as.list(solution_terms(ar_search))[0:suggested_size])
  )
  if(!silent){message(paste0("Done!\nP_perp = ", P_perp))}
  
  # extract residuals from projected AR component
  projected_ar <- projpred::project(ar_search, nterms = suggested_size)
  proj_ar_pred <- colMeans(projpred::proj_linpred(projected_ar)$pred)
  res <- proj_ar_pred - head(data, length(proj_ar_pred))
  
  # apply projpred to MA component
  if(!silent){message("Applying projpred to the MA component ... ")}
  ma_refmodel <- .fit_ma_model(res, Q = Q_max)
  ma_search <- .perform_search(ma_refmodel)
  Q_perp <- projpred::suggest_size(
    ma_search, alpha = alpha, pct = pct, stat = stat
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
