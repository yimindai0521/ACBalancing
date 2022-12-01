#' @title Univariate Balancing
#' @description Univariate Balancing is a univariate perspective of
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
#' @inheritParams MB
#'
#' \code{delta.space} grid of values for the tuning parameter, a vector of
#' candidate values for the degree of approximate covariate balance. The default
#' is c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6).
#'
#'
#' @return a UB object with the following attributes:
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
#'
#' @rdname UB
#' @export
#'
UB <- function(covariate, treat, group1, group2 = NULL, outcome, opti.method = "proximal",
               delta.space = c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
               iterations = 1000, convergence = 10^{-10}, rate = 0.0001, bootstrap.time = 1000) {

  # Initiazing input
  covariate <- as.matrix(covariate)
  treat <- as.vector(treat)
  outcome <- as.vector(outcome)
  delta.space <- as.vector(delta.space)

  # Check the compatibility of covariate, treat and outcome
  try(if (nrow(covariate) != length(treat)) stop("the sample size of covariate and treat are not the same!"))
  try(if (nrow(covariate) != length(outcome)) stop("the sample size of covariate and treat are not the same!"))
  try(if (length(treat) != length(outcome)) stop("the sample size of covariate and treat are not the same!"))

  # Check the compatibility of group1/group2 and treat
  try(if (length(group1) > 1) stop("the length of group1 should be one!"))
  try(if (sum(treat == group1) == 0) stop("the sample size of treat = group1 is zero!"))
  if(!is.null(group2)){
    try(if (sum(treat == group2) == 0) stop("the sample size of treat = group2 is zero!"))
  }

  # Initiaze delta.space and Check the compatibility of delta.space
  delta.space <- unique(delta.space)
  try(if (sum(delta.space >= 0) == 0) stop("all the value in delta.space is negative!"))

  # Check the compatibility of method and opti.method
  stopifnot(opti.method %in% c("proximal", "BFGS"))

  # Check the compatibility of iterations, convergence and rate
  try(if (iterations <= 0) stop("iterations should be greater than zero!"))
  try(if (convergence < 0) stop("convergence should be non-negative!"))
  try(if (rate < 0) stop("rate should be non-negative!"))

  # Initiazing information for data
  covariate <- scale(covariate)
  data.matrix <- as.matrix(cbind(covariate, treat, outcome))
  dimension <- dim(covariate)[2]
  sample.size <- dim(covariate)[1]
  group1.number <- sum(treat == group1)
  xMB <- CholMB(x = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB", dimension = dimension)

  # Set up space for saving weight and GMIM
  weight.space <- matrix(0, nrow = group1.number, ncol = length(delta.space))
  GMIM <- matrix(0, nrow = bootstrap.time, ncol = length(delta.space))
  x <- t(covariate[treat == group1, ])
  # Calculate weight for each delta in delta.space
  for (i in 1:length(delta.space)) {
    eps <- 1000
    beta <- as.vector(rep(0, dimension))

    # Implementation of BFGS method
    if (opti.method == "BFGS") {
      Opti.func <- function(beta) {
        sum(exp(t(x) %*% beta - 1)) + delta.space[i] * norm(beta, type = c("1"))
      }
      beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))$par
    }

    # Implementation of proximal optimization method with R code
    if (opti.method == "proximal") {
      iterationtime = 0;
      while (eps > convergence) {
        beta_old <- beta
        hessian <- solve(crossprod(t(x), t(x) * matrix(exp(crossprod(x, beta_old)), ncol(x), nrow(x))) + rate * diag(rep(1, nrow(x))))
        vt <- beta - crossprod(t(hessian), crossprod(t(x), exp(crossprod(x, beta))))
        beta <- soft(vt, delta.space[i], norm = "l1")
        fobj_old <- sum(exp(t(x) %*% beta_old)) + delta.space[i] * sum(abs(beta_old))
        fobj_new <- sum(exp(t(x) %*% beta)) + delta.space[i] * sum(abs(beta))
        eps <- abs(fobj_old - fobj_new)
        if(iterationtime >= 1000){
          eps = 0
        }
      }
    }

    weight <- exp(crossprod(x, beta))
    weight <- weight / sum(weight)
    weight.space[, i] <- weight

    # Bootstrap to select tuning parameter
    for(j in 1:bootstrap.time){
      index <- sample(1:group1.number, replace = TRUE)
      GMIM[j, i] <- sum(tcrossprod(t(weight[index]), x[, index])^2)
    }
  }

  # calculate the AT
  meanGMIM <- colMeans(GMIM)
  delta <- delta.space[rank(meanGMIM) == 1]
  weight <- weight.space[, rank(meanGMIM) == 1]
  GMIM <- sum(tcrossprod(t(weight), xMB)^2)
  AT <- crossprod(weight, outcome[treat == group1])
  result <- list(AT = AT, weight = weight, GMIM = GMIM, delta = delta, parameter = -beta)
  return(result)
}
