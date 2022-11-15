#' @title Mahalanobis balancing
#' @description Mahalanobis balancing is a multivariate perspective of
#' approximate covariate balancing method to estimate average treatment
#' effect.
#'
#' @details
#' \code{group1} and \code{group0}
#' \itemize{
#' \item To estimate \eqn{E(Y (1))} (average treatment effect for group 1),
#' you need to set \code{group1} = 1 and ignore \code{group2}. Similarly,
#' To estimate \eqn{E(Y (0))} (average treatment effect for group 0),
#' you need to set \code{group1} = 0 and ignore \code{group2}.
#'
#' \item To estimate average treatment effect on the control group \eqn{E(Y (1) | T = 0)},
#' you need to set \code{group1} = 1 and \code{group2} = 0. Similarly, To estimate
#' average treatment effect on the treated group \eqn{E(Y (0) | T = 1)},
#' you need to set \code{group1} = 0 and \code{group2} = 1.
#'
#' \item This function is feasible when there are more than two groups.
#' }
#'
#' \code{method} can be a valid string, including
#' \itemize{
#' \item "MB": We choose the weighting matrix \eqn{{W}_1=[diag(\hat{\Sigma})]^{-1}}
#' where \eqn{\hat{\Sigma}} denotes sample covariance matrix.
#' \item "MB2": We choose the weighting matrix \eqn{{W}_2={[\hat{\Sigma}]^{-1}}}
#' where \eqn{\hat{\Sigma}} denotes sample covariance matrix.
#' \item "kernelMB": Firstly, we modify our covariate to \eqn{X_i* = (\Phi(X_1,X_i), ..., \Phi(X_n,X_i))},
#' then we apply method "MB" to produce Mahalanobis balancing weights.
#' }
#'
#' \code{delta.space} grid of values for the tuning parameter, a vector of
#' candidate values for the degree of approximate covariate balance. The default
#' is c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6).
#'
#' @param x covariates.
#' @param treat treatment indicator vector.
#' @param group1 see Details.
#' @param group2 see Details, Default: NA.
#' @param outcome outcome vector.
#' @param method a string that takes values in {"MB", "MB2", "kernelMB"}. See Details. Default: 'MB'.
#' @param delta.space tuning parameter in balancing. See Details. Default: c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6).
#' @param iterations iteration time in optimization problem, Default: 1000.
#'
#' @return a MB object with the following attributes:
#' \itemize{
#' \item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
#' \item{weight:}{ the estimated Mahalanobis balancing weight.}
#' \item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper.}
#' \item{delta:}{ the tuning parameter we choose.}
#' }
#'
#' @examples
#' ## estimating ATE##
#' set.seed(0521)
#' data <- si.data()
#' result1 <- MB(x = data$X, treat = data$Tr, group1 = 1, outcome = data$Y, method = "MB")
#' result2 <- MB(x = data$X, treat = data$Tr, group1 = 0, outcome = data$Y, method = "MB")
#'
#' ## an estimate of ATE
#' result1$AT - result2$AT
#'
#' ## estimating ATC##
#' result3 <- MB(x = data$X, treat = data$Tr, group1 = 1, group2 = 0, outcome = data$Y, method = "MB")
#'
#' ## an estimate of ATC
#' result3$AT - mean(data$Y[data$Tr == 0])
#'
#' @rdname MB
#' @export
#'
MB <- function(x, treat, group1, group2 = NA, outcome, method = "MB", delta.space = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), iterations = 1000) {
  # data
  stopifnot(method %in% c("MB", "MB2", "kernelMB"))

  data.matrix <- cbind(x, treat, outcome)
  dimension <- dim(x)[2]
  sample.size <- dim(x)[1]
  if (is.null(dimension)) {
    dimension <- 1
    sample.size <- sum(abs(x) >= 0)
  }
  group1.number <- sum(treat == group1)

  if (method == "MB") {
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = group2, method = "MB", dimension = dimension)
  }
  if (method == "MB2") {
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = group2, method = "MB2", dimension = dimension)
  }
  if (method == "kernelMB") {
    dimension <- dim(x)[1]
    data <- matrix(NA, dim(x)[1], dim(x)[1])
    data.temp <- matrix(NA, dim(x)[1], dim(x)[1])
    for (i in 1:sample.size) {
      for (j in 1:sample.size) {
        data.temp[i, j] <- sum((data.matrix[i, 1:dim(x)[2]] - data.matrix[j, 1:dim(x)[2]])^2)
      }
    }
    median.scale <- median(data.temp)
    for (i in 1:sample.size) {
      for (j in 1:sample.size) {
        data[i, j] <- exp(-data.temp[i, j] / median.scale)
      }
    }
    solution <- Cholesky(x = x, treat = treat, group1 = group1, group2 = group2, method = "MB2", dimension = dimension)
  }
  x <- solution

  # calculate the weights for group1
  w <- rep(0, group1.number)
  u1 <- rep(0, dimension)
  GMIM <- rep(NA, sum(delta.space > 0))
  for (i in 1:sum(delta.space > 0)) {
    Opti.func <- function(u) {
      u1 <- t(x) %*% u
      sum(exp(u1 - 1)) + sqrt(delta.space[i]) * norm(u, type = c("2"))
    }
    solution <- optim(par = u, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))
    u.sol <- t(x) %*% solution$par
    w <- exp(u.sol - 1)
    w <- w / sum(w)
    GMIM[i] <- t(w) %*% t(x) %*% x %*% w
  }

  # calculate the AT
  delta <- delta.space[rank(GMIM) == 1]
  Opti.func <- function(u) {
    u1 <- t(x) %*% u
    sum(exp(u1 - 1)) + sqrt(delta) * norm(u, type = c("2"))
  }
  solution <- optim(par = u, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))
  u.sol <- t(x) %*% solution$par
  w <- exp(u.sol - 1)
  w <- w / sum(w)
  GMIM <- t(w) %*% t(x) %*% x %*% w

  dimension <- dim(x)[1]
  AT <- t(w) %*% outcome[treat == group1]
  result <- list(AT = AT, weight = w, GMIM = GMIM, delta = delta, parameter = -solution$par)
  return(result)
}
