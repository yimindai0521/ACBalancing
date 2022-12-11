#' @title Univariate Balancing
#' @description Univariate Balancing (Wang, Y., & Zubizarreta, J. R. (2020)) is a univariate perspective of
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
#' \code{delta.space} grid of values for the tuning parameter, a vector of
#' candidate values for the degree of approximate univariate covariate balance. The default
#' is c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1).
#'
#' @inheritParams MB
#' @param opti.method The optimization method to optimize loss function. It should takes values in {"BFGS", "proximal"}. Default: "proximal".
#' @param bootstrap.time The bootstrap time to select tuning parameter.
#'
#' @return a UB object with the following attributes:
#' \itemize{
#' \item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
#' \item{weight:}{ the estimated univariate balancing weight.}
#' \item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper.}
#' \item{delta:}{ the tuning parameter we choose.}
#' }
#'
#' @rdname UB
#'
#' @references
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately
#' balancing weights: asymptotic properties and practical considerations,
#' \emph{Biometrika}
#' \strong{107(1), 93-105.},
#' \doi{10.1093/biomet/asz050}.
#'
#' @export
#'
#' @examples
#' # estimating ATE
#' set.seed(0521)
#' data <- si.data()
#' X <- data$X
#' treat <- data$Tr
#' Y <- data$Y
#'
#' result1 <- UB(covariate = X, treat = treat, group1 = 1, outcome = Y)
#' result2 <- UB(covariate = X, treat = treat, group1 = 0, outcome = Y)
#'
#' # an estimate of ATE
#' result1$AT - result2$AT
#'
#' # estimating ATC
#' result3 <- UB(covariate = X, treat = treat, group1 = 1, group2 = 0, outcome = Y)
#'
#' # an estimate of ATC
#' result3$AT - mean(data$Y[data$Tr == 0])
#'
UB <- function(covariate, treat, group1, group2 = NULL, outcome, opti.method = c("BFGS", "proximal"),
               delta.space = c(1e-04, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1),
               iterations = 1000, convergence = 1e-8, rate = 10, bootstrap.time = 1000) {

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
  if(!is.null(group2)){
    try(if (sum(treat == group2) == 0) stop("the sample size of treat = group2 is zero!"))
  }

  # initialize delta.space and Check the compatibility of delta.space
  delta.space <- unique(delta.space)
  try(if (sum(delta.space >= 0) == 0) stop("all the value in delta.space is negative!"))

  # Check the compatibility of method and opti.method
  stopifnot(opti.method %in% c("proximal", "BFGS"))

  # Check the compatibility of iterations, convergence and rate
  try(if (iterations <= 0) stop("iterations should be greater than zero!"))
  try(if (convergence < 0) stop("convergence should be non-negative!"))
  try(if (rate < 0) stop("rate should be non-negative!"))

  # initialize information for data
  covariate <- scale(covariate) # Standardize covariate
  dimension <- dim(covariate)[2] # Calculate the dimension of covariate
  sample.size <- dim(covariate)[1] # Calculate the sample size of covariate
  group1.number <- sum(treat == group1) # Calculate the number of individuals in group1
  xMB <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB") # Calculate the scaled covariate using diagonal matrix

  # Set up space for saving weight and GMIM
  weight.space <- matrix(0, nrow = group1.number, ncol = length(delta.space))
  beta.space <- matrix(0, nrow = dimension, ncol = length(delta.space))
  GMIM <- matrix(0, nrow = bootstrap.time, ncol = length(delta.space))

  # Select the individual in group1
  x <- t(covariate[treat == group1, ])

  # Calculate weight for each delta in delta.space
  for (i in 1:length(delta.space)) {
    eps <- 1000 # initiazing eps
    beta <- as.vector(rep(0, dimension)) # initiazing beta

    opti.method <- match.arg(opti.method) # Choose optimization method
    # Implementation of BFGS method
    if (opti.method == "BFGS") {
      Opti.func <- function(beta) {
        sum(exp(t(x) %*% beta - 1)) + delta.space[i] * sum(abs(beta)) # loss function for univariate balancing
      }
      beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))$par # Using BFGS to optimize loss function
    }

    # Implementation of proximal optimization method with R code
    if (opti.method == "proximal") {
      iterationtime = 0 # initialize iteration time
      while (eps > convergence) {
        beta_old <- beta # Save beta before iteration
        hessian <- solve(crossprod(t(x), t(x) * matrix(exp(crossprod(x, beta_old)), ncol(x), nrow(x))) + rate * diag(rep(1, nrow(x)))) # Calculate modified Hessian Newton update with parameter lambda > 0
        vt <- beta - crossprod(t(hessian), crossprod(t(x), exp(crossprod(x, beta)))) # Gradient step
        beta <- soft(vt, delta.space[i], norm = "l1") # Proximal operator step
        fobj_old <- sum(exp(t(x) %*% beta_old)) + delta.space[i] * sum(abs(beta_old)) # Calculate the value of loss function with beta_old
        fobj_new <- sum(exp(t(x) %*% beta)) + delta.space[i] * sum(abs(beta)) # Calculate the value of loss function with beta
        eps <- abs(fobj_old - fobj_new)  # Calculate the difference of the value
        if(iterationtime >= iterations){
          eps = 0 # If the iterationtime is greater /equal to iterations, then stop optimization
        }
      }
    }

    # Calculate weight
    weight <- exp(crossprod(x, beta)) # Calculate weight
    weight <- weight / sum(weight) # Normalize weight

    # Save result in corresponding space
    beta.space[, i] <- beta
    weight.space[, i] <- weight

    # Bootstrap to select tuning parameter
    for(j in 1:bootstrap.time){
      index <- sample(1:group1.number, replace = TRUE) # sampling index
      GMIM[j, i] <- sum(tcrossprod(t(weight[index]), x[, index])^2) # Calculate GMIM in bootstraped sample
    }
  }

  # calculate the AT
  meanGMIM <- colMeans(GMIM) # Calculate the average error in bootstrap.time
  delta <- delta.space[which.min(meanGMIM)] # Select the tuning parameter in loss function
  beta <- as.vector(beta.space[, which.min(meanGMIM)]) # Find the solution to parameter with selected tuning parameter
  weight <- weight.space[, which.min(meanGMIM)] # Find the solution to weight with selected tuning parameter
  GMIM <- sum(tcrossprod(t(weight), xMB)^2) # Calculate GMIM
  AT <- crossprod(weight, outcome[treat == group1]) # Calculate average treatment effect in target group
  result <- list(AT = AT, weight = as.vector(weight), GMIM = GMIM, delta = delta, parameter = as.vector(-beta))
  return(result)
}
