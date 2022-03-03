

#' Construction of distance covariance optimal weights weights
#'
#' @description Constructs independence-inducing weights (distance covariance optimal weights) for 
#' estimation of causal quantities for continuous-valued treatments
#'
#' @param A vector indicating the value of the treatment or exposure variable. Should be a numeric vector.
#' @param X matrix of covariates with number of rows equal to the length of \code{A} and each column is a
#' \strong{pre-treatment} covariate to be balanced between treatment groups.
#' @param lambda tuning parameter for the penalty on the sum of squares of the weights
#' @param decorrelate_moments logical scalar. Whether or not to add constraints that result in exact decorrelation of 
#' weighted first order moments of \code{X} and \code{A}. Defaults to \code{FALSE}.
#' @param preserve_means logical scalar. Whether or not to add constraints that result in exact preservation of 
#' weighted first order moments of \code{X} and \code{A}. Defaults to \code{FALSE}.
#' @param dimension_adj logical scalar. Whether or not to add adjustment to energy distance terms that account for 
#' the dimensionality of \code{X}. Defaults to \code{TRUE}.
#' @return An object of class \code{"independence_weights"} with elements:
#' \item{weights}{A vector of length \code{nrow(X)} containing the estimated sample weights }
#' \item{A}{Treatment vector}
#' \item{opt}{The optimization object returned by \code{osqp::solve_osqp()}}
#' \item{objective}{The value of the objective function at its optimal value. This is the weighted dependence statistic plus any ridge penalty on the weights.}
#' \item{D_unweighted}{The value of the weighted dependence distance using all weights = 1 (i.e. unweighted)}
#' \item{D_w}{The value of the weighted dependence distance of Huling, et al. (2021) using the optimal estimated weights. This is the weighted dependence statistic without the ridge penalty on the weights.}
#' \item{distcov_unweighted}{The unweighted distance covariance term. This is the standard distance covariance of Szekely et al (2007). This term
#' is always equal to \code{D_unweighted}.}
#' \item{distcov_weighted}{The weighted distance covariance term. This term itself does not directly measure weighted dependence but is a critical component of it.  }
#' \item{energy_A}{The weighted energy distance between \code{A} and its weighted version}
#' \item{energy_X}{The weighted energy distance between \code{X} and its weighted version}
#' \item{ess}{The estimated effective sample size of the weights using Kish's effective sample size formula.}
#' @seealso \code{\link[independenceWeights]{print.independence_weights}} for printing of fitted energy balancing objects
#' @importFrom osqp solve_osqp
#' @references Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. 
#' Annals of Statistics 35(6) 2769-2794 \url{https://doi.org/10.1214/009053607000000505}
#' 
#' Huling, J. D., Greifer, N., & Chen, G. (2021). Independence weights for causal inference with continuous exposures. 
#' arXiv preprint arXiv:2107.07086. \url{https://arxiv.org/abs/2107.07086}
#'
#' @examples
#'
#' simdat <- simulate_confounded_data(seed = 999, nobs = 500)
#'
#' y <- simdat$data$Y
#' A <- simdat$data$A
#' X <- as.matrix(simdat$data[c("Z1", "Z2", "Z3", "Z4", "Z5")])
#'
#' dcows <- independence_weights(A, X)
#'
#' print(dcows)
#'
#' # distribution of response:
#' quantile(y)
#' 
#' ## create grid
#' trt_vec <- seq(min(simdat$data$A), 50, length.out=500)
#' 
#' ## estimate ADRF
#' adrf_hat <- weighted_kernel_est(A, y, dcows$weights, trt_vec)$est
#' 
#' ## estimate naively without weights
#' adrf_hat_unwtd <- weighted_kernel_est(A, y, rep(1, length(y)), trt_vec)$est
#' 
#' ylims <- range(c(simdat$data$Y, simdat$true_adrf(trt_vec)))
#' plot(x = simdat$data$A, y = simdat$data$Y, ylim = ylims, xlim = c(0,50))
#' ## true ADRF
#' lines(x = trt_vec, y = simdat$true_adrf(trt_vec), col = "blue", lwd=2)
#' ## estimated ADRF
#' lines(x = trt_vec, y = adrf_hat, col = "red", lwd=2)
#' ## naive estimate
#' lines(x = trt_vec, y = adrf_hat_unwtd, col = "green", lwd=2)
#'
#' @export
independence_weights <- function(A, 
                                 X, 
                                 lambda = 0, 
                                 decorrelate_moments = FALSE,
                                 preserve_means = FALSE,
                                 dimension_adj = TRUE)
{
  
  weights <- rep(1, NROW(A))
  
  n <- NROW(A)
  p <- NCOL(X)
  gamma <- 1
  
  stopifnot(gamma >= 0)
  
  Xdist  <- as.matrix(dist(X))
  Adist  <- as.matrix(dist(A))
  
  stopifnot(n == NROW(X))
  
  ## terms for energy-dist(Wtd A, A)
  Q_energy_A  <- -Adist / n ^ 2
  aa_energy_A <- 1 * as.vector(rowSums(Adist)) / (n ^ 2)
  
  ## terms for energy-dist(Wtd X, X)
  Q_energy_X  <- -Xdist / n ^ 2
  aa_energy_X <- 1 * as.vector(rowSums(Xdist)) / (n ^ 2)
  
  
  mean_Adist <- mean(Adist)
  mean_Xdist <- mean(Xdist)
  
  
  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")
  
  
  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")
  
  ## quadratic term for weighted total distance covariance
  P <- XA * AA / n ^ 2
  
  
  #Constraints: positive weights, weights sum to n
  if (preserve_means)
  {
    
    if (decorrelate_moments)
    {
      Constr_mat <- drop(scale(A, scale=FALSE)) * scale(X, scale=FALSE)
      Amat <- rbind(diag(n), rep(1, n), t(X), A, t(Constr_mat))
      
      lvec <- c(rep(0, n), n, colMeans(X), mean(A), rep(0, ncol(X)))
      uvec <- c(rep(Inf, n), n, colMeans(X), mean(A), rep(0, ncol(X)))
    } else
    {
      Amat <- rbind(diag(n), rep(1, n), t(X), A)
      
      lvec <- c(rep(0, n), n, colMeans(X), mean(A))
      uvec <- c(rep(Inf, n), n, colMeans(X), mean(A))
    }
    
  } else 
  {
    if (decorrelate_moments)
    {
      #Constr_mat <- (A - mean(A)) * scale(X, scale = FALSE)
      Constr_mat <- drop(scale(A, scale=FALSE)) * scale(X, scale=FALSE)
      
      Amat <- rbind(diag(n), rep(1, n), t(Constr_mat))
      lvec <- c(rep(0, n), n, rep(0, ncol(X)))
      uvec <- c(rep(Inf, n), n, rep(0, ncol(X)))
    } else
    {
      Amat <- rbind(diag(n), rep(1, n))
      lvec <- c(rep(0, n), n)
      uvec <- c(rep(Inf, n), n)
    }
  }
  
  
  if (dimension_adj)
  {
    Q_energy_A_adj <- 1 / sqrt(p)
    Q_energy_X_adj <- 1
    
    sum_adj <- 1*(Q_energy_A_adj + Q_energy_X_adj)
    
    Q_energy_A_adj <- Q_energy_A_adj / sum_adj
    Q_energy_X_adj <- Q_energy_X_adj / sum_adj
    
  } else
  {
    Q_energy_A_adj <- Q_energy_X_adj <- 1
  }
  
  #Optimize. try up to 15 times until there isn't a weird failure of solve_osqp()
  for (na in 1:15)
  {
    opt.out <- solve_osqp(2 * (P + gamma * (Q_energy_A * Q_energy_A_adj + Q_energy_X * Q_energy_X_adj) + lambda * diag(n) / n ^ 2 ),
                          q = 2 * gamma * (aa_energy_A * Q_energy_A_adj + aa_energy_X * Q_energy_X_adj),
                          A = Amat, l = lvec, u = uvec,
                          pars = osqp::osqpSettings(max_iter = 2e5,
                                                    eps_abs = 1e-8,
                                                    eps_rel = 1e-8,
                                                    verbose = FALSE))
    
    if (!identical(opt.out$info$status, "maximum iterations reached") & !(any(opt.out$x > 1e5)))
    {
      break
    }
  }
  
  weights <- opt.out$x
  weights[weights < 0] <- 0 #due to numerical imprecision
  
  ## quadratic part of the overall objective function
  QM_unpen <- (P + gamma * (Q_energy_A * Q_energy_A_adj + Q_energy_X * Q_energy_X_adj) )
  #QM <- QM_unpen + lambda * diag(n)
  quadpart_unpen <- drop(t(weights) %*% QM_unpen %*% weights)
  quadpart_unweighted <- sum(QM_unpen)
  quadpart <- quadpart_unpen + sum(weights ^ 2) * lambda / n ^ 2
  
  ## linear part of the overall objective function
  qvec <- 2 * gamma * (aa_energy_A * Q_energy_A_adj + aa_energy_X * Q_energy_X_adj)
  linpart  <- drop(weights %*% qvec)
  linpart_unweighted  <- sum(qvec)
  
  ## objective function
  objective_history <- quadpart + linpart + gamma*(-1*mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj)
  
  
  ## D(w)
  D_w <- quadpart_unpen + linpart + gamma*(-1*mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj)
  
  D_unweighted <- quadpart_unweighted + linpart_unweighted + gamma*(-1*mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj)
  
  
  qvec_full <- 2 * (aa_energy_A * Q_energy_A_adj + aa_energy_X * Q_energy_X_adj)
  
  #quadpart_energy <- drop(t(weights) %*% ((Q_energy_A + Q_energy_X)) %*% weights)
  
  quadpart_energy_A <- drop(t(weights) %*% ((Q_energy_A)) %*% weights) * Q_energy_A_adj
  quadpart_energy_X <- drop(t(weights) %*% ((Q_energy_X)) %*% weights) * Q_energy_X_adj
  
  quadpart_energy <- quadpart_energy_A * Q_energy_A_adj + quadpart_energy_X * Q_energy_X_adj
  
  
  distcov_history <- drop(t(weights) %*% P %*% weights)
  
  unweighted_dist_cov <- sum(P)
  
  
  linpart_energy   <- drop(weights %*% qvec_full)
  linpart_energy_A <- 2 * drop(weights %*% (aa_energy_A)) * Q_energy_A_adj
  linpart_energy_X <- 2 * drop(weights %*% (aa_energy_X)) * Q_energy_X_adj
  
  
  ## sum of energy-dist(wtd A, A)+energy-dist(wtd X, X)
  energy_history <- quadpart_energy + linpart_energy - mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj
  
  ## energy-dist(wtd A, A)
  energy_A <- quadpart_energy_A + linpart_energy_A - mean_Adist * Q_energy_A_adj
  
  ## energy-dist(wtd X, X)
  energy_X <- quadpart_energy_X + linpart_energy_X - mean_Xdist * Q_energy_X_adj
  
  
  objective_history <- objective_history
  energy_history    <- energy_history
  
  ess <- (sum(weights)) ^ 2 / sum(weights ^ 2)
  
  ret_obj <- list(weights = weights, 
                  A = A,
                  opt = opt.out, 
                  objective = objective_history,     ### the actual objective function value
                  D_unweighted = D_unweighted,
                  D_w = D_w,
                  distcov_unweighted = unweighted_dist_cov,
                  distcov_weighted = distcov_history,         ### the weighted total distance covariance
                  energy_A = energy_A,               ### Energy(Wtd Treatment, Treatment)
                  energy_X = energy_X,               ### Energy(Wtd X, X)
                  ess = ess)        # effective sample size
  
  class(ret_obj) <- c("independence_weights")
  
  return(ret_obj)              
}



#' Calculation of weighted energy statistics for weighted dependence
#'
#' @description Calculates weighted energy statistics used to quantify weighted dependence
#'
#' @param A treatment vector indicating values of the treatment/exposure variable.
#' @param X matrix of covariates with number of rows equal to the length of \code{weights} and each column is a
#' covariate
#' @param weights a vector of sample weights
#' @param dimension_adj logical scalar. Whether or not to add adjustment to energy distance terms that account for 
#' the dimensionality of \code{x}. Defaults to \code{TRUE}.
#' @return a list with the following components
#' \item{D_w}{The value of the weighted dependence distance of Huling, et al. (2021) using the optimal estimated weights. 
#' This is the weighted dependence statistic without the ridge penalty on the weights.}
#' \item{distcov_unweighted}{The unweighted distance covariance term. This is the standard distance covariance of Szekely et al (2007). This term
#' is always equal to \code{D_unweighted}.}
#' \item{distcov_weighted}{The weighted distance covariance term. This term itself does not directly measure weighted dependence but is a critical component of it.  }
#' \item{energy_A}{The weighted energy distance between \code{A} and its weighted version}
#' \item{energy_X}{The weighted energy distance between \code{X} and its weighted version}
#' \item{ess}{The estimated effective sample size of the weights using Kish's effective sample size formula.}
#' @references Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances. 
#' Annals of Statistics 35(6) 2769-2794 \url{https://doi.org/10.1214/009053607000000505}
#' 
#' Huling, J. D., Greifer, N., & Chen, G. (2021). Independence weights for causal inference with continuous exposures. 
#' arXiv preprint arXiv:2107.07086. \url{https://arxiv.org/abs/2107.07086}
#' @examples 
#' 
#' simdat <- simulate_confounded_data(seed = 999, nobs = 100)
#' 
#' str(simdat$data)
#' 
#' A <- simdat$data$A
#' X <- as.matrix(simdat$data[c("Z1", "Z2", "Z3", "Z4", "Z5")])
#' 
#' wts <- runif(length(A))
#' 
#' weighted_energy_stats(A, X, wts)
#'
#' @export
weighted_energy_stats <- function(A, X, weights,
                                  dimension_adj = TRUE)
{
  
  Xdist  <- as.matrix(dist(X))
  Adist  <- as.matrix(dist(A))
  
  gamma <- 1
  
  ## normalize weights
  weights <- weights / mean(weights)
  
  n <- NROW(A)
  p <- NCOL(X)
  
  ## terms for energy-dist(Wtd A, A)
  Q_energy_A  <- -Adist / n ^ 2
  aa_energy_A <- 1 * as.vector(rowSums(Adist)) / (n ^ 2)
  
  ## terms for energy-dist(Wtd X, X)
  Q_energy_X  <- -Xdist / n ^ 2
  aa_energy_X <- 1 * as.vector(rowSums(Xdist)) / (n ^ 2)
  
  
  mean_Adist <- mean(Adist)
  mean_Xdist <- mean(Xdist)
  
  
  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")
  
  
  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")
  
  ## quadratic term for weighted total distance covariance
  P <- XA * AA / n ^ 2
  
  
  
  
  if (dimension_adj)
  {
    Q_energy_A_adj <- 1 / sqrt(p)
    Q_energy_X_adj <- 1
    
    sum_adj <- 1*(Q_energy_A_adj + Q_energy_X_adj)
    
    Q_energy_A_adj <- Q_energy_A_adj / sum_adj
    Q_energy_X_adj <- Q_energy_X_adj / sum_adj
    
  } else
  {
    Q_energy_A_adj <- Q_energy_X_adj <- 1
  }
  
  
  ## quadratic part of the overall objective function
  QM <- (P + gamma * (Q_energy_A * Q_energy_A_adj + Q_energy_X * Q_energy_X_adj) )
  quadpart <- drop(t(weights) %*% QM %*% weights)
  
  ## linear part of the overall objective function
  qvec <- 2 * gamma * (aa_energy_A * Q_energy_A_adj + aa_energy_X * Q_energy_X_adj)
  linpart  <- drop(weights %*% qvec)
  
  ## objective function
  objective_history <- quadpart + linpart + gamma*(-1*mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj)
  
  
  qvec_full <- 2 * (aa_energy_A * Q_energy_A_adj + aa_energy_X * Q_energy_X_adj)
  
  
  quadpart_energy_A <- drop(t(weights) %*% ((Q_energy_A)) %*% weights) * Q_energy_A_adj
  quadpart_energy_X <- drop(t(weights) %*% ((Q_energy_X)) %*% weights) * Q_energy_X_adj
  
  quadpart_energy <- quadpart_energy_A * Q_energy_A_adj + quadpart_energy_X * Q_energy_X_adj
  
  
  distcov_history <- drop(t(weights) %*% P %*% weights)
  
  
  
  linpart_energy   <- drop(weights %*% qvec_full)
  linpart_energy_A <- 2 * drop(weights %*% (aa_energy_A)) * Q_energy_A_adj
  linpart_energy_X <- 2 * drop(weights %*% (aa_energy_X)) * Q_energy_X_adj
  
  
  ## sum of energy-dist(wtd A, A)+energy-dist(wtd X, X)
  energy_history <- quadpart_energy + linpart_energy - mean_Xdist * Q_energy_X_adj - mean_Adist * Q_energy_A_adj
  
  ## energy-dist(wtd A, A)
  energy_A <- quadpart_energy_A + linpart_energy_A - mean_Adist * Q_energy_A_adj
  
  ## energy-dist(wtd X, X)
  energy_X <- quadpart_energy_X + linpart_energy_X - mean_Xdist * Q_energy_X_adj

  objective_history <- objective_history
  energy_history    <- energy_history
  
  ess <- (sum(weights)) ^ 2 / sum(weights ^ 2)
  
  retobj <- list(D_w = objective_history,           ### the actual objective function value
                 distcov_unweighted = sum(P),
                 distcov_weighted = distcov_history,         ### the weighted total distance covariance
                 energy_A = energy_A,               ### Energy(Wtd Treatment, Treatment)
                 energy_X = energy_X,               ### Energy(Wtd X, X)
                 ess = ess)                         ### effective sample size
  
  class(retobj) <- c("weighted_energy_terms", "list")
  
  return(retobj)
}

