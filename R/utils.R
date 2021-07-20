

#' Calculation of weighted nonparametric regression estimate of the dose response function
#'
#' @description Calculates weighted nonparametric regression estimate of the causal average dose response function
#'
#' @param A vector indicating the value of the treatment or exposure variable. Should be a numeric vector.
#' @param y vector of responses
#' @param weights a vector of sample weights of length equal to the length of \code{y}
#' @param Aseq a vector of new points for which to obtain estimates of E(Y(a))
#' @return a vector of estimates of a causal ADRF at the values of the treatment specified by \code{Aseq}
#'
#' @export
weighted_kernel_est <- function(A, y, weights, Aseq)
{
  
  stopifnot(length(A) == length(y))
  stopifnot(length(A) == length(weights))
  
  ## if the weights are degenerate, do nothing!
  if (all(weights == 0) | (mean(weights == 0) > 0.95))
  {
    estimated <- rep(NA, NROW(Aseq))
    locpoly_fit <- NULL
  } else {
    dfx <- data.frame(Y = y, TRT = A)
    locpoly_fit <- locfit(y ~ lp(A), weights = weights, data = dfx)
    
    estimated <- unname(locfit:::predict.locfit(locpoly_fit, 
                                                newdata = Aseq,
                                                where = "fitp"))
  }
  
  estimated
}
  