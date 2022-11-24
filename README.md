# ACBalancing

An R package for Approximate Covariate Balancing.

<!-- badges: start -->
<!-- badges: end -->

## Overview

The R package ACBalancing implements several approximate balancing 
  methods, such as the univariate approximate balancing method in [Minimal dispersion approximately balancing weights: asymptotic properties and practical considerations](https://doi.org/10.1093/biomet/asz050), the Mahalanobis
  balancing method in [Mahalanobis balancing: a multivariate perspective on approximate covariate balancing](https://arxiv.org/abs/2204.13439), and the hierarchically regularized entropy balancing method in [Hierarchically Regularized Entropy Balancing](http://dx.doi.org/10.2139/ssrn.3807620),
  for estimating average treatment effects (ATE) or average treatment effects in 
  controlled or treated groups (ATC or ATT). It includes Mahalanobis balancing and high-dimensional Mahalanobis balancing method. 


## Installation

You can install the development version of ACBalancing from [GitHub](https://github.com/) with:

``` r
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("yimindai0521/ACBalancing")
```

## Usage Examples

We illustrate the usage of MBalance package using simple synthetic datasets.

``` r

```

## To do list

Implement function UB, HRB (including add reference for function), add comments on code, move the optimization code to C++ and complete vignettes.
