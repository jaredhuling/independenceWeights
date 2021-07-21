


#' Simulation of confounded data with a continuous treatment
#'
#' @description Simulates confounded data with continuous treatment based on Vegetabile et al's simulation
#'
#' @param seed random seed for reproducibility
#' @param nobs number of observations
#' @param MX1 the mean of the first covariate. Defaults to -0.5, the value used in the simulations of Vegetabile, et al (2021).
#' @param MX2 the mean of the second and fourth covariates. Defaults to 1, the value used in the simulations of Vegetabile, et al (2021).
#' @param MX3 the probability that the fifth covariate (a binary covariate) is equal to 1. Defaults to 0.3, the value used in the 
#' simulations of Vegetabile, et al (2021).
#' @param A_effect whether (\code{TRUE}) or not (\code{FALSE}) the treatment has a causal effect on the outcome. If \code{TRUE}, the 
#' setting used is that of the main text of Vegetabile, et al (2021). If \code{FALSE}, the setting is that used in the Appendix of
#' Vegetabile, et al (2021).
#' @return An list with elements:
#' \item{data}{A simulated dataset with \code{nobs} rows }
#' \item{true_adrf}{A function that inputs values of the treatment \code{A} and outputs the true ADRF, E(Y(A)), of the data-generating
#' mechanism used to generate \code{data}. }
#' @references Vegetabile, B. G., Griffin, B. A., Coffman, D. L., Cefalu, M., Robbins, M. W., and McCaffrey, D. F. (2021). 
#' Nonparametric estimation of population average dose-response curves using entropy balancing weights for continuous exposures. 
#' Health Services and Outcomes Research Methodology, 21(1), 69-110.
#' @importFrom stats dist rbinom rchisq rnorm
#'
#' @examples
#' 
#' simdat <- simulate_confounded_data(seed = 999, nobs = 500)
#' 
#' str(simdat$data)
#' 
#' A <- simdat$data$A
#' y <- simdat$data$Y
#' 
#' trt_vec <- seq(min(simdat$data$A), max(simdat$data$A), length.out=500)
#' ylims <- range(c(simdat$data$Y, simdat$true_adrf(trt_vec)))
#' plot(x = simdat$data$A, y = simdat$data$Y, ylim = ylims)
#' lines(x = trt_vec, y = simdat$true_adrf(trt_vec), col = "blue", lwd=2)
#' 
#' ## naive estimate of ADRF without weights
#' adrf_hat_unwtd <- weighted_kernel_est(A, y, rep(1, length(y)), trt_vec)
#' lines(x = trt_vec, y = adrf_hat_unwtd, col = "green", lwd=2)
#' 
#' 
#' @export
simulate_confounded_data <- function(seed = 1,
                                     nobs = 1000,
                                     MX1 =  -0.5,
                                     MX2  = 1,
                                     MX3 =  0.3,
                                     A_effect = TRUE)
{
  
  set.seed(seed)
  
  ## setup parameters
  # MX1 =  -0.5
  #MX2  = 1
  # MX3 =  0.3
  # a_effect =TRUE
  ## generate covariates and dose
  
  X1 <- rnorm(nobs, mean = MX1, sd = 1)
  X2 <- rnorm(nobs, mean = MX2, sd = 1)
  X3 <- rnorm(nobs, mean = 0, sd = 1)
  X4 <- rnorm(nobs, mean = MX2, sd = 1)
  X5 <- rbinom(nobs, 1, prob = MX3)
  
  Z1 <- exp(X1 / 2)
  Z2 <- (X2 / (1 + exp(X1))) + 10
  Z3 <- (X1 * X3 / 25) + 0.6
  Z4 <- (X4 - MX2)**2 
  Z5 <- X5
  
  
  muA <- 5 * abs(X1) + 6 * abs(X2) + 3 * abs(X5) + abs(X4)
  
  A <- rchisq(nobs, df = 3, ncp = muA)
  
  
  if(A_effect)
  {
    Cnum <- ((MX1+3)^2+1) + 2*((MX2-25)^2+1)
    Y <- - 0.15 * A^2 + A * (X1^2 + X2^2) - 15 + (X1+3)^2 + 2 * (X2-25)^2 + X3 - Cnum + rnorm(nobs, sd = 1)
    Y <- Y / 50
    truth <- - 0.15 * A^2 + A * (2 + MX1^2 + MX2^2) - 15 #+ 5 * (1 + MX1^2 + 6 * MX1 + 9) + 15 * (1 + MX2^2 + 6 * MX2 + 9) + MX3
    truth <- truth / 50
    # truth <- - (A_test - 10) * (A_test - 10) / 5 + 5* A_test * (2 + MX1^2 + MX2^2) / 1 - 15 * (MX1 - MX2) + MX3 + (1 + MX2^2 - 40 * MX2 + 400)
  } else 
  {
    Y <- X1 + X1^2 + X2 + X2^2 + X1 * X2 + X5 + rnorm(nobs, sd = 1)
    truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3, nobs)
  }
  
  ## true treatment effect curve
  Afunc <- function(A)
  {
    if(A_effect)
    {
      truth <- - 0.15 * A^2 + A * (2 + MX1^2 + MX2^2) - 15
      truth <- truth / 50
    } else 
    {
      truth <- rep(MX1 + (MX1^2 + 1) + MX2 + (MX2^2 + 1) + MX1 * MX2 + MX3, length(A))
    }
    return(truth)
  }
  
  datz <- data.frame('Y' = Y, 'A' = A, 'Z1' = Z1, 'Z2' = Z2, 'Z3' = Z3, 'Z4' = Z4, 'Z5' = Z5, 'truth' = truth)
  
  datx <- data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5)
  
  list(data = datz, true_adrf = Afunc, original_covariates = datx)
}




