#' Hierarchically Regularized Entropy Balancing
#' @description Hierarchically Regularized Entropy Balancing is a
#' approximate balancing method that expands the feature space by including
#' higher-order terms of covariates while imposes
#' ridge penalties with a hierarchical structure on the higher-order terms.
#'
#' @inheritParams MB
#'
#' @param alpha.space alpha.space is
#' @param cv.fold cv.fold is
#' @param second.moment second.moment is
#' @param third.moment third.moment is
#' @param interact interact is
#'
#' @return a HRB object with the following attributes:
#' \itemize{
#' \item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
#' \item{weight:}{ the estimated univariate balancing weight.}
#' \item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper.}
#' \item{MAE:}{ Mean Absolute Error.}
#' \item{alpha:}{ the tuning parameter we choose.}
#' }
#'
#' @rdname HRB
#'
#' @export
#'
#' @examples
#' # estimating ATE
#' set.seed(0521)
#' data <- si.data()
#' result1 <- HRB(x = data$X, treat = data$Tr, group1 = 1, outcome = data$Y)
#' result2 <- HRB(x = data$X, treat = data$Tr, group1 = 0, outcome = data$Y)
#'
#' # an estimate of ATE
#' result1$AT - result2$AT
#'
#' # estimating ATC
#' result3 <- MB(x = data$X, treat = data$Tr, group1 = 1, group2 = 0, outcome = data$Y, method = "MB")
#'
#' # an estimate of ATC
#' result3$AT - mean(data$Y[data$Tr == 0])
HRB <- function(covariate, treat, group1, group2 = NULL, outcome,
                alpha.space = c(1e-1, 0.05, 1e-2, 0.005, 1e-3),
                iterations = 1000, convergence = 10^{-10}, cv.fold = 5,
                second.moment = TRUE, third.moment = TRUE, interact = TRUE) {

  # Initiazing input
  covariate <- as.matrix(covariate)
  treat <- as.vector(treat)
  outcome <- as.vector(outcome)
  alpha.space <- as.vector(alpha.space)

  # Adding higher-order terms of covariates
  new_covariate <- covariate
  if (second.moment == TRUE) {
    new_covariate <- cbind(new_covariate, covariate^2)
  }
  if (interact == TRUE) {
    for(k in 1:ncol(covariate)){
      new_covariate <- cbind(new_covariate, covariate[, k] * covariate[, (1:ncol(covariate) > k)])
    }
  }
  if (third.moment == TRUE) {
    new_covariate <- cbind(new_covariate, covariate^3)
  }
  covariate <- new_covariate

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

  # Initiaze alpha.space and Check the compatibility of alpha.space
  alpha.space <- unique(alpha.space)
  try(if (sum(alpha.space >= 0) == 0) stop("all the value in alpha.spaceis negative!"))

  # Check the compatibility of iterations, convergence and rate
  try(if (iterations <= 0) stop("iterations should be greater than zero!"))
  try(if (convergence < 0) stop("convergence should be non-negative!"))

  # Initiazing information for data
  covariate <- scale(covariate)
  data.matrix <- as.matrix(cbind(covariate, treat, outcome))
  dimension <- dim(covariate)[2]
  sample.size <- dim(covariate)[1]
  group1.number <- sum(treat == group1)
  xMB <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB")

  # Set up space for saving weight and GMIM
  weight.space <- matrix(0, nrow = group1.number, ncol = length(alpha.space))
  beta.space <- matrix(0, nrow = dimension, ncol = length(alpha.space))
  GMIM <- matrix(0, nrow = cv.fold, ncol = length(alpha.space))
  x <- t(covariate[treat == group1, ])

  #
  if (is.null(group2)) {
    M <- colMeans(covariate)
  } else {
    M <- colMeans(covariate[treat == group2, ])
  }

  # Calculate weight for each delta in alpha.space
  fold_ids <- sample(1:group1.number) %% cv.fold + 1
  MAE <- matrix(NA, length(alpha.space), cv.fold)
  for (i in 1:length(alpha.space)) {
    for (j in 1:cv.fold) {
      eps <- 1000
      beta <- as.vector(rep(0, dimension))
      x.train <- x[, fold_ids != j]
      x.test <- x[, fold_ids == j]
      Opti.func <- function(beta) {
        log(sum(exp(-crossprod(x.train, beta)))) + crossprod(M, beta) + alpha.space[i] * norm(beta, type = c("2"))^2
      }
      beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = convergence, maxit = iterations))$par

      # Bootstrap to select tuning parameter
      weight.test <- exp(-crossprod(x.test, beta))
      weight.test <- weight.test / sum(weight.test)

      #
      MAE[j, i] <- mean(abs(tcrossprod(t(weight.test), x.test)))
    }
  }

  # calculate the AT
  meanMAE <- colMeans(MAE)
  alpha <- alpha.space[which.min(meanMAE)]

  eps <- 1000
  beta <- as.vector(rep(0, dimension))
  Opti.func <- function(beta) {
    log(sum(exp(-crossprod(x, beta)))) + crossprod(M, beta) + alpha * norm(beta, type = c("2"))^2
  }
  beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))$par
  weight <- exp(-crossprod(x, beta))
  weight <- weight / sum(weight)
  GMIM <- sum(tcrossprod(t(weight), xMB)^2)
  AT <- crossprod(weight, outcome[treat == group1])
  result <- list(AT = AT, weight = as.vector(weight), GMIM = GMIM, alpha = alpha, parameter = as.vector(beta))
  return(result)
}
