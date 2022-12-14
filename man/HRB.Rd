% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRB.R
\name{HRB}
\alias{HRB}
\title{Hierarchically Regularized Entropy Balancing}
\usage{
HRB(
  covariate,
  treat,
  group1,
  group2 = NULL,
  outcome,
  alpha.space = c(0.1, 0.05, 0.01, 0.005, 0.001),
  iterations = 1000,
  convergence = 1e-08,
  cv.fold = 5,
  second.moment = TRUE,
  third.moment = TRUE,
  interact = TRUE
)
}
\arguments{
\item{covariate}{Covariates, it should be a n by d matrix where n is sample size and d is dimension.}

\item{treat}{The vector of treatment assignments, it should be a vector of length n.}

\item{group1}{group1 is a number, see Details.}

\item{group2}{group2 is a number, see Details, Default: NULL.}

\item{outcome}{The outcome vector, it should be a vector of length n.}

\item{alpha.space}{alpha.space is grid of values for the tuning parameter, a vector of
candidate values for the degree of the strength of penalty.}

\item{iterations}{The maximum number of iteration times in optimization problem, Default: 1000.}

\item{convergence}{The absolute tolerance used to determine a stopping criteria for the optimization problems. Default: 10^{-8}.}

\item{cv.fold}{cv.fold is the fold of cross validation. It should be a integer.}

\item{second.moment}{second.moment is a bool variable to judge whether you will add second moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.}

\item{third.moment}{third.moment is a bool variable to judge whether you will add third moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.}

\item{interact}{interact is a bool variable to judge whether you will add interaction moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.}
}
\value{
a HRB object with the following attributes:
\itemize{
\item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
\item{weight:}{ the estimated univariate balancing weight.}
\item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper (Dai, Y., & Yan, Y. (2022).).}
\item{MAE:}{ Mean Absolute Error that defines in (Xu, Y., & Yang, E. (2021)).}
\item{alpha:}{ the tuning parameter we choose.}
}
}
\description{
Hierarchically Regularized Entropy Balancing is a
approximate balancing method that expands the feature space by including
higher-order terms of covariates while imposes
ridge penalties with a hierarchical structure on the higher-order terms.
}
\examples{
# estimating ATE
set.seed(0521)
data <- si.data()
X <- data$X
treat <- data$Tr
Y <- data$Y

result1 <- HRB(covariate = X, treat = treat, group1 = 1, outcome = Y)
result2 <- HRB(covariate = X, treat = treat, group1 = 0, outcome = Y)

# an estimate of ATE
result1$AT - result2$AT

# estimating ATC
result3 <- HRB(covariate = X, treat = treat, group1 = 1, group2 = 0, outcome = Y)

# an estimate of ATC
result3$AT - mean(data$Y[data$Tr == 0])
}
\references{
Xu, Y., & Yang, E. (2021). Hierarchically Regularized Entropy Balancing. \emph{Political Analysis}, \strong{forthcoming}, \doi{10.2139/ssrn.3807620}.

Dai, Y., & Yan, Y. (2022). Mahalanobis balancing: a multivariate perspective on approximate covariate balancing. arXiv preprint arXiv:2204.13439.
}
