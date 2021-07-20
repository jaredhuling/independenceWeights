
#' Printing results for estimated energy balancing weights
#'
#' @description Prints results for energy balancing weights
#'
#' @param x a fitted object from \code{\link[independenceWeights]{independence_weights}}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @rdname print
#' @seealso \code{\link[independenceWeights]{independence_weights}} for function which produces energy balancing weights
#' @importFrom stats quantile
#' @export
print.independence_weights <- function(x, digits = max(getOption('digits')-3, 3), ...)
{
  
  cat(#"Unweighted dependence distance: ", round(x$D_unweighted, digits),
      "Unweighted distance covariance:          ", round(x$distcov_unweighted, digits),
      "\nOptimized weighted dependence distance:  ", round(x$D_w, digits),
      "\nEffective sample size:                 ", round(x$ess, digits), "\n\n")
  
  cat("Weight ranges:\n")
  print(summary(x$weights, digits = digits), digits = digits)
}




#' Printing weighted energy statistics for given weights
#'
#' @description Prints weighted energy statistics for given weights
#'
#' @param x a fitted object from \code{\link[independenceWeights]{weighted_energy_stats}}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @rdname print
#' @seealso \code{\link[independenceWeights]{weighted_energy_stats}} for function which produces energy balancing weights
#' @importFrom stats quantile
#' @export
print.weighted_energy_terms <- function(x, digits = max(getOption('digits')-3, 3), ...)
{
  
  cat(#"Unweighted dependence distance: ", round(x$D_unweighted, digits),
    "Unweighted distance covariance:           ", round(x$distcov_unweighted, digits),
    "\nWeighted dependence distance:             ", round(x$D_w, digits),
    "\nWeighted energy distance(A, weighted A):  ", round(x$energy_A, digits),
    "\nWeighted energy distance(X, weighted X):  ", round(x$energy_X, digits),
    "\nEffective sample size:                  ", round(x$ess, digits), "\n\n")
}
