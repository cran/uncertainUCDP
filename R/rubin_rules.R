#' Pool regression coefficients and standard errors via Rubin's rules
#'
#' Pools coefficient estimates and standard errors from multiple fitted models
#' (typically one per dataset) using Rubin's rules, returning a
#' `coeftest`-like matrix with pooled estimates, SEs, t-statistics, and p-values. Can be used when when pooling mulitple models with different draws from the UCDP uncertainty model.
#'
#' If `cluster` is provided, cluster-robust covariance matrices are computed for
#' each model via `sandwich::vcovCL()` and passed to `lmtest::coeftest()`.
#'
#' @param list A list of fitted model objects (e.g., `lm`, `glm`, etc.)
#'   for which `lmtest::coeftest()` methods exist.
#' @param cluster Optional. A clustering variable passed to `sandwich::vcovCL()`. Typically a one-sided formula specficying clusters, e.g., `~ cluster_id`
#' @param df Optional. Degrees of freedom for t-distribution when computing p-values. If `NULL`, uses the degrees of freedom from the first model's `coeftest` output, i.e. the df from the original model.
#'
#' @return A numeric matrix with the same layout as `lmtest::coeftest()`,
#'   containing pooled coefficients, pooled standard errors, t-statistics, and
#'   two-sided p-values. The returned object uses the first element's
#'   `coeftest` matrix as a template (including row/column names and attributes).
#'
#' @details
#' The pooling follows the usual within/between decomposition:
#' within-imputation variance is the mean of squared SEs, and between-imputation
#' variance is the sample variance of point estimates across imputations. The
#' total variance is \eqn{W + (1 + 1/m)B} where \eqn{m} is the number of models.
#'
#' @examples
#' library(foreach)
#' library(dplyr)
#' # Draw 100 datasets from UCDP uncertainty model and return as a list
#' data('ucdpged')
#' ucdp_sb <- ucdpged %>% filter(type_of_violence == 1)
#' list_of_models <- foreach(i = 1:100) %do%{
#'  lm(ged_sb_draw ~ year + region, data =
#' ucdp_sb %>% mutate(ged_sb_draw = runcertainUCDP(n = n(), fatalities = best, tov = 'sb')))
#' }
#' #' # Pool model coefficients and standard errors via Rubin's rules
#' pooled_results <- models_rubin_rules(list_of_models)
#'
#'
#'
#' @export
models_rubin_rules <- function(list, df = NULL, cluster = NULL) {
  if (!is.null(cluster)) {
    clusters <- foreach::foreach(i = 1:length(list)) %do%
      {
        sandwich::vcovCL(list[[i]], cluster = cluster)
      }
    coeftests <- foreach::foreach(i = 1:length(list)) %do%
      {
        lmtest::coeftest(list[[i]], vcov = clusters[[i]], df = df)
      }
  } else {
    coeftests <- foreach::foreach(i = 1:length(list)) %do%
      {
        lmtest::coeftest(list[[i]], df = df)
      }
  }

  out <- rubin_b_and_ses(coeftests)

  return(out)
}

#' Pool coefficients and standard errors from multiple `coeftest` objects
#'
#' Internal helper that takes a list of `lmtest::coeftest()` matrices and returns
#' a pooled `coeftest`-like matrix using Rubin's rules.
#'
#' @param coeftests A list of objects returned by `lmtest::coeftest()`. Each
#'   element must have matching coefficient rows/columns, where column 1 is the
#'   estimate and column 2 is the standard error.
#'
#' @return A `coeftest`-like numeric matrix with pooled estimates (col 1), pooled
#'   standard errors (col 2), pooled t-statistics (col 3), and two-sided p-values
#'   (col 4). Row/column names and attributes are inherited from the first
#'   element of `coeftests`.
#'
rubin_b_and_ses <- function(coeftests) {
  coefs <- foreach::foreach(
    i = 1:length(coeftests),
    .final = dplyr::bind_rows
  ) %do%
    {
      coeftests[[i]][, 1]
    }

  ses <- foreach::foreach(
    i = 1:length(coeftests),
    .final = dplyr::bind_rows
  ) %do%
    {
      coeftests[[i]][, 2]
    }

  coef_means <- colMeans(coefs)
  coef_sqdiff <- foreach::foreach(i = 1:ncol(coefs), .final = unlist) %do%
    {
      sum((coefs[, i] - coef_means[i])^2)
    }

  ses_sq <- foreach::foreach(i = 1:ncol(ses), .final = unlist) %do%
    {
      sum(ses[, i]^2)
    }

  ses <- sqrt(
    ses_sq / nrow(ses) + (1 + 1 / nrow(ses)) * coef_sqdiff / (length(coefs) - 1)
  )

  ts <- coef_means / ses
  pvals <- 2 * stats::pt(-abs(ts), df = attributes(coeftests[[1]])$df)

  out <- coeftests[[1]]
  out[, 1] <- coef_means
  out[, 2] <- ses
  out[, 3] <- ts
  out[, 4] <- pvals
  return(out)
}
