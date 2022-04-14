sarma.sim = function(n, phi1, phi2, theta1, theta2){
  n = 100
  y = numeric(n)
  y[1] = y[2] = 0
  eps = rnorm(n)
  
  for(i in 3:n)
  {
    y[i] = phi1*y[i-1] + phi2*y[i-2] + theta1*eps[i-1] + theta2*eps[i-2] + eps[i]
  }
  
  return(y)
}

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
    ar_fit <- brms::brm(
        ar_formula,
        data = X,
        family = gaussian,
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
    ma_fit <- brms::brm(
        ma_formula,
        data = E,
        family = gaussian,
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

.get_ma_refmodel <- function(data, P_perp, Q_max = 6){
    # fit restricted AR model in BRMS and extract residuals
    ar_fit <- .fit_ar_model(data, P_perp)
    res <- residuals(ar_fit)
    eps <- res[,1]
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
        vs <- projpred::varsel(refmodel, method = "arma", verbose = FALSE)
    }
    else {
        vs <- projpred::cv_varsel(refmodel, method = "arma", verbose = FALSE)
    }
    # use heuristic to chose informative parameters in submodel
    suggested_size <- suggest_size(vs, alpha = alpha, pct = pct, stat = stat)
    # extract the largest lag from the restricted covariate space
    lag_perp <- length(
        unlist(as.list(projpred::solution_terms(vs))[0:suggested_size])
    )
    # lag_perp <- length(cov_perp) %>%
    #    strsplit(".t") %>%
    #    sapply("[", 2) %>%
    #    strtoi() %>%
    #    max()
    return(lag_perp)
}

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
    P_perp <- .get_restricted_lag(
        ar_refmodel, data, alpha = alpha, pct = pct, cv = cv, stat = stat
    )
    if(!silent){message(paste0("Done!\nP_perp = ", P_perp))}
    # apply projpred to MA component
    if(!silent){message("Applying projpred to the MA component ... ")}
    ma_refmodel <- .get_ma_refmodel(data, P_perp = P_perp, Q_max = Q_max)
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
