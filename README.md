# ACBalancing

An R package for Approximate Covariate Balancing.

<!-- badges: start -->

<!-- badges: end -->

## Overview

The R package ACBalancing implements several approximate balancing methods, such as the univariate approximate balancing method in [Minimal dispersion approximately balancing weights: asymptotic properties and practical considerations](https://doi.org/10.1093/biomet/asz050), the Mahalanobis balancing method in [Mahalanobis balancing: a multivariate perspective on approximate covariate balancing](https://arxiv.org/abs/2204.13439), and the hierarchically regularized entropy balancing method in [Hierarchically Regularized Entropy Balancing](http://dx.doi.org/10.2139/ssrn.3807620), for estimating average treatment effects (ATE) or average treatment effects in controlled or treated groups (ATC or ATT). It includes Mahalanobis balancing and high-dimensional Mahalanobis balancing method.

## Installation

You can install the development version of ACBalancing from [GitHub](https://github.com/) with:

``` r
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("yimindai0521/ACBalancing")
```

## Usage Examples

We illustrate the usage of MBalance package using simple datasets.

``` r
set.seed(0521)
data <- si.data()
result1 <- MB(covariate = data$X, treat = data$Tr, group1 = 1, outcome = data$Y, opti.method = "proximalC", method = "MB")
result2 <- MB(covariate = data$X, treat = data$Tr, group1 = 0, outcome = data$Y, opti.method = "proximalC", method = "MB")
result3 <- UB(covariate = data$X, treat = data$Tr, group1 = 1, outcome = data$Y, opti.method = "proximal")
result4 <- UB(covariate = data$X, treat = data$Tr, group1 = 0, outcome = data$Y, opti.method = "proximal")
result5 <- HRB(covariate = data$X, treat = data$Tr, group1 = 1, outcome = data$Y, second.moment = FALSE, third.moment = FALSE, interact = FALSE)
result6 <- HRB(covariate = data$X, treat = data$Tr, group1 = 0, outcome = data$Y, second.moment = FALSE, third.moment = FALSE, interact = FALSE)

# estimating ATE
result1$AT - result2$AT # produce an estimate of ATE using the Mahalanobis balancing method
result3$AT - result4$AT # produce an estimate of ATE using the Univariate Balancing method
result5$AT - result6$AT # produce an estimate of ATE using the Hierarchically Regularized Entropy Balancing method

# parameter for propensity score (\beta_1-\beta_0)
result1$parameter - result2$parameter # produce an estimate of parameter in propensity score using the Mahalanobis balancing method
result3$parameter - result4$parameter # produce an estimate of parameter in propensity score using the Univariate Balancing method
result5$parameter - result6$parameter # produce an estimate of parameter in propensity score using the Hierarchically Regularized Entropy Balancing method

# Generalized Mahalobnis Imbalance Measure
result1$GMIM + result2$GMIM # Generalized Mahalanobis Imbalance Measure for Mahalanobis balancing weights
result3$GMIM + result4$GMIM # Generalized Mahalanobis Imbalance Measure for Univariate Balancing weights
result5$GMIM + result6$GMIM # Generalized Mahalanobis Imbalance Measure for Hierarchically Regularized Entropy Balancing weights
```

For more details, please check website: https://yimindai0521.github.io/ACBalancing.
