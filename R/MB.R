#' @title Mahalanobis balancing
#' @description Mahalanobis balancing (Dai, Y., & Yan, Y. (2022).) is a multivariate perspective of
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
#' candidate values for the degree of multivariate approximate covariate balance.
#' The default is c(1e-1, 0.05, 1e-2, 0.005, 1e-3).
#'
#' @param covariate Covariates, it should be a n by d matrix where n is sample size and d is dimension.
#' @param treat The vector of treatment assignments, it should be a vector of length n.
#' @param group1 group1 is a number, see Details.
#' @param group2 group2 is a number, see Details, Default: NULL.
#' @param outcome The outcome vector, it should be a vector of length n.
#' @param opti.method The optimization method that will be applied to optimize loss function. It should takes values in {"BFGS", "proximal", "proximalC"}. Default: "proximalC".
#' @param method Method to calculate weighting matrix. It should takes values in {"MB", "MB2", "kernelMB"}. See Details. Default: "MB".
#' @param delta.space Tuning parameter in balancing. See Details. Default: c(1e-1, 0.05, 1e-2, 0.005, 1e-3).
#' @param convergence The absolute tolerance used to determine a stopping criteria for the optimization problems. Default: 10^{-8}.
#' @param rate The parameter in Tikhonov-regularized Newton update. Default: 10.
#' @param iterations The maximum number of iteration times in optimization problem, Default: 1000.
#'
#' @return a MB object with the following attributes:
#' \itemize{
#' \item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
#' \item{weight:}{ the estimated Mahalanobis balancing weight.}
#' \item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper.}
#' \item{delta:}{ the tuning parameter we choose.}
#' }
#'
#' @references Dai, Y., & Yan, Y. (2022). Mahalanobis balancing: a multivariate perspective on approximate covariate balancing. arXiv preprint arXiv:2204.13439.
#'
#' @examples
#' # estimating ATE
#' set.seed(0521)
#' data <- si.data()
#' X <- data$X
#' treat <- data$Tr
#' Y <- data$Y
#'
#' result1 <- MB(covariate = X, treat = treat, group1 = 1, outcome = Y)
#' result2 <- MB(covariate = X, treat = treat, group1 = 0, outcome = Y)
#'
#' # an estimate of ATE
#' result1$AT - result2$AT
#'
#' # estimating ATC
#' result3 <- MB(covariate = X, treat = treat, group1 = 1, group2 = 0, outcome = Y)
#'
#' # an estimate of ATC
#' result3$AT - mean(data$Y[data$Tr == 0])
#'
#' @rdname MB
#' @importFrom stats optim cov dist median rnorm runif
#' @export
MB <- function(covariate, treat, group1, group2 = NULL, outcome, opti.method = c("proximalC", "proximal", "BFGS"),
               method = c("MB", "MB2", "kernelMB"), delta.space = c(1e-1, 0.05, 1e-2, 0.005, 1e-3), iterations = 1000,
               convergence = 1e-8, rate = 10) {

  # initialize input
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
  if (!is.null(group2)) {
    try(if (sum(treat == group2) == 0) stop("the sample size of treat = group2 is zero!"))
  }

  # Initiaze delta.space and Check the compatibility of delta.space
  delta.space <- unique(delta.space)
  try(if (sum(delta.space >= 0) == 0) stop("all the value in delta.space is negative!"))

  # Check the compatibility of method and opti.method
  stopifnot(method %in% c("MB", "MB2", "kernelMB"))
  stopifnot(opti.method %in% c("proximalC", "proximal", "BFGS"))

  # Check the compatibility of iterations, convergence and rate
  try(if (iterations <= 0) stop("iterations should be greater than zero!"))
  try(if (convergence < 0) stop("convergence should be non-negative!"))
  try(if (rate < 0) stop("rate should be non-negative!"))

  # Initiazing information for data
  data.matrix <- as.matrix(cbind(covariate, treat, outcome)) # Save the covariate and treat in the same matrix
  dimension <- dim(covariate)[2] # Calculate the dimension of covariate
  sample.size <- dim(covariate)[1] # Calculate the sample size of covariate
  group1.number <- sum(treat == group1) # Calculate the number of covariate

  # Calculate the scaled covariate
  method <- match.arg(method)
  if (method == "MB") {
    x <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB") # Calculate the scaled covariate using diagonal matrix
  }
  if (method == "MB2") {
    x <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB2") # Calculate the scaled covariate using covariance matrix
  }
  if (method == "kernelMB") {
    dimension <- dim(covariate)[1] # re-calculate the dimension
    data <- matrix(NA, dim(covariate)[1], dim(covariate)[1]) # Open a space to save the new covariate
    dist.matrix <- as.matrix(dist(covariate)) # Open a space to save distance matrix
    median.scale <- median(dist.matrix) # Calculate the bandwidth for kernel function
    covariate <- exp(-dist.matrix / median.scale) # Calculate the transformed covariate
    x <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB") # Calculate the scaled covariate after transfromation using diagonal matrix
  }

  # Set up space for saving weight, beta and GMIM
  GMIM <- rep(NA, sum(delta.space > 0))
  beta.space <- matrix(0, nrow = dimension, ncol = length(delta.space))
  weight.space <- matrix(0, nrow = group1.number, ncol = length(delta.space))

  # Calculate weight for each delta in delta.space
  for (i in 1:length(delta.space)) {
    eps <- 1000 # initialize eps
    beta <- as.vector(rep(0, dimension)) # initialize beta

    # Implementation of BFGS method
    opti.method <- match.arg(opti.method)
    if (opti.method == "BFGS") {
      Opti.func <- function(beta) {
        sum(exp(t(x) %*% beta - 1)) + delta.space[i] * norm(beta, type = c("2")) # Loss function to be optimized
      }
      beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))$par # Using BFGS to optimize loss function
    }

    # Implementation of proximal optimization method with R code
    if (opti.method == "proximal") {
      iterationtime <- 0 # Initiazing iteration time
      while (eps > convergence) {
        beta_old <- beta # Save beta before iteration
        hessian <- solve(crossprod(t(x), t(x) * matrix(exp(crossprod(x, beta_old)), ncol(x), nrow(x))) + rate * diag(rep(1, nrow(x)))) # Calculate modified Hessian Newton update with parameter lambda > 0
        vt <- beta - crossprod(t(hessian), crossprod(t(x), exp(crossprod(x, beta)))) # Gradient step
        beta <- soft(vt, delta.space[i]) # Proximal operator step
        fobj_old <- sum(exp(t(x) %*% beta_old)) + delta.space[i] * sqrt(sum(beta_old^2)) # Calculate the value of loss function with beta_old
        fobj_new <- sum(exp(t(x) %*% beta)) + delta.space[i] * sqrt(sum(beta^2)) # Calculate the value of loss function with beta
        eps <- abs(fobj_old - fobj_new) # Calculate the difference of the value
        iterationtime <- iterationtime + 1
        if (iterationtime >= iterations) {
          eps <- 0 # If the iterationtime is greater /equal to iterations, then stop optimization
        }
      }
    }

    # Implementation of proximal optimization method with C++ code
    if (opti.method == "proximalC") {
      beta <- proximaloptim(X = x, rate = rate, lambda = delta.space[i], beta_start = beta, convergence = convergence, iteration = iterations)
    }

    # Calculate weight and balancing diagnostics
    w <- exp(crossprod(x, beta)) # Calculate weight
    w <- w / sum(w) # Normalize weight
    GMIM[i] <- sum(tcrossprod(t(w), x)^2) # Calculate Generalized Mahalanobis Imbalance Measure

    # Save weight and beta to corresponding space
    beta.space[, i] <- beta
    weight.space[, i] <- w
  }

  # calculate the AT
  delta <- as.vector(delta.space[which.min(GMIM)]) # Select the tuning parameter in loss function
  GMIM <- min(GMIM) # Find the GMIM with selected tuning parameter
  beta <- beta.space[, which.min(GMIM)] # Find the solution to parameter with selected tuning parameter
  weight <- weight.space[, which.min(GMIM)] # Find the solution to weight with selected tuning parameter
  AT <- crossprod(weight, outcome[treat == group1]) # Calculate average treatment effect in target group
  result <- list(AT = AT, weight = as.vector(weight), GMIM = GMIM, delta = delta, parameter = as.vector(-beta)) # list the result and return
  return(result)
}
