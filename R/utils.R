

#' Calculation of weighted nonparametric regression estimate of the dose response function
#'
#' @description Calculates weighted nonparametric regression estimate of the causal average dose response function
#'
#' @param A vector indicating the value of the treatment or exposure variable. Should be a numeric vector.
#' @param y vector of responses
#' @param weights a vector of sample weights of length equal to the length of \code{y}
#' @param Aseq a vector of new points for which to obtain estimates of E(Y(a))
#' @return A list with the following elements
#' \item{fit}{A fitted model object from the \code{\link[locfit]{lp}} function}
#' \item{estimated}{a vector of estimates of a causal ADRF at the values of the treatment specified by \code{Aseq}}
#' @importFrom stats predict
#'
#' @export
weighted_kernel_est <- function(A, y, weights, Aseq)
{
  
  stopifnot(length(A) == length(y))
  if (missing(weights)) weights <- rep(1, length(A))
  stopifnot(length(A) == length(weights))
  
  if (missing(Aseq)) {
    min_A <- min(A)
    max_A <- max(A)
    range_A <- max_A - min_A
    Aseq <- seq(min_A - .05*range_A, max_A + .05*range_A, length.out = 500)
  }
  
  ## if the weights are degenerate, do nothing!
  if (all(weights == 0) || (mean(weights == 0) > 0.95))
  {
    estimated <- rep(NA_real_, length(Aseq))
    locpoly_fit <- NULL
  } else {
    dfx <- data.frame(Y = y, TRT = A)
    locpoly_fit <- locfit::locfit(Y ~ locfit::lp(TRT), weights = weights, data = dfx)
    
    estimated <- unname(predict(locpoly_fit, 
                                newdata = Aseq,
                                where = "fitp"))
  }
  
  list(fit = locpoly_fit, estimated = estimated)
}
  
