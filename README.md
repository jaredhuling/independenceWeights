
<!-- README.md is generated from README.Rmd. Please edit that file -->

# independenceWeights

<!-- badges: start -->
<!-- badges: end -->

The goal of independenceWeights is to …

## Installation

You can install the released version of independenceWeights from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("independenceWeights")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jaredhuling/independenceWeights")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(independenceWeights)
#> Loading required package: osqp
## basic example code
```

Simulate data with a continuous treatment that has a confounded
relationship with a response:

``` r
simdat <- simulate_confounded_data(seed = 999, nobs = 500)
y <- simdat$data$Y ## response
A <- simdat$data$A ## treatment
X <- as.matrix(simdat$data[c("Z1", "Z2", "Z3", "Z4", "Z5")]) ## confounders
```

Now estimate weights to adjust for confounders using the distance
covariance optimal weights (DCOWs), which aim to mitigate the dependence
between *A* and *X*:

``` r
dcows <- independence_weights(A, X)

print(dcows)
#> Unweighted distance covariance:          0.3963 
#> Optimized weighted dependence distance:  0.0246 
#> 
#> Weight ranges:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  0.2767  0.8215  1.0000  1.4120  5.7360
```

Now use the weights to estimate the causal average dose response
function (ADRF)

``` r
## create grid
trt_vec <- seq(min(simdat$data$A), 50, length.out=500)

## estimate ADRF
adrf_hat <- weighted_kernel_est(A, y, dcows$weights, trt_vec)

## estimate naively without weights
adrf_hat_unwtd <- weighted_kernel_est(A, y, rep(1, length(y)), trt_vec)

ylims <- range(c(simdat$data$Y, simdat$true_adrf(trt_vec)))
plot(x = simdat$data$A, y = simdat$data$Y, ylim = ylims, xlim = c(0,50),
     xlab = "A", ylab = "Y")
## true ADRF
lines(x = trt_vec, y = simdat$true_adrf(trt_vec), col = "blue", lwd=2)
## estimated ADRF
lines(x = trt_vec, y = adrf_hat, col = "red", lwd=2)
## naive estimate
lines(x = trt_vec, y = adrf_hat_unwtd, col = "green", lwd=2)
```

<img src="man/figures/README-adrf-1.png" width="100%" />
