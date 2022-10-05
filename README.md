# Kullback-Leibler Projections for Automatic Order Identification of Seasonal ARMA Models

This repository holds the code for the paper _Bayesian order identification of ARMA models with projection predictive inference_ by McLatchie et al. ([2022](https://arxiv.org/abs/2208.14824)).

## Algorithm

The primary algorithm presented in the paper is as follows:

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

This code uses a local branch of `projpred` in order to perform the time series search heuristic, where the only change from vanilla `projpred` is in `/R/search.R` where we add
```diff
search_forward <- function(p_ref, refmodel, nterms_max, verbose = TRUE, opt,
               submodls = submodels))
 }


+search_arma <- function(p_ref, refmodel, nterms_max, verbose = TRUE, opt,
+                           search_terms = NULL) {
+  iq <- ceiling(quantile(seq_len(nterms_max), 1:10 / 10))
+  if (is.null(search_terms)) {
+    allterms <- split_formula(refmodel$formula, data = refmodel$fetch_data())
+  } else {
+    allterms <- search_terms
+  }
+
+  chosen <- character()
+  total_terms <- count_terms_chosen(allterms)
+  stop_search <- min(total_terms, nterms_max)
+  submodels <- c()
+
+  for (size in seq_len(stop_search)) {
+    cands <- select_possible_terms_size(chosen, allterms, size = size)
+    # select the next lag
+    cands <- cands[1]
+    if (is.null(cands))
+      next
+    full_cands <- lapply(cands, function(cand) c(chosen, cand))
+    subL <- lapply(full_cands, project_submodel,
+                   p_ref = p_ref, refmodel = refmodel, regul = opt$regul)
+
+    ## select best candidate
+    imin <- which.min(sapply(subL, "[[", "kl"))
+    chosen <- c(chosen, cands[imin])
+
+    ## append submodels
+    submodels <- c(submodels, list(subL[[imin]]$submodl))
+
+    if (verbose && count_terms_chosen(chosen) %in% iq) {
+      print(paste0(names(iq)[max(which(count_terms_chosen(chosen) == iq))],
+                   " of terms selected."))
+    }
+  }
+
+  ## reduce chosen to a list of non-redundant accumulated models
+  return(list(solution_terms = setdiff(reduce_models(chosen), "1"),
+              submodls = submodels))
+}
+
```

in `R/varsel.R` where we add

```diff
select <- function(method, p_sel, refmodel, nterms_max, penalty, verbose, opt,
                                   search_terms = search_terms)
     search_path$p_sel <- p_sel
     return(search_path)

+} else if (method == "arma") {
+    search_path <- search_arma(p_sel, refmodel, nterms_max, verbose, opt,
+                               search_terms = search_terms)
+    search_path$p_sel <- p_sel
+    return(search_path)
}
}

...

parse_args_varsel <- function(refmodel, method, refit_prj, nterms_max,
     }
   }

-  if (!(method %in% c("l1", "forward"))) {
+  if (!(method %in% c("l1", "forward", "arma"))) {
     stop("Unknown search method")
   }

```
