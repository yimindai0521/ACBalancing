#' Hierarchically Regularized Entropy Balancing
#' @description Hierarchically Regularized Entropy Balancing is a
#' approximate balancing method that expands the feature space by including
#' higher-order terms of covariates while imposes
#' ridge penalties with a hierarchical structure on the higher-order terms.
#'
#' @details
#' \code{group1} and \code{group0}
#' \itemize{
#' \item group1 and group0 are not continuous.
#'
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
#' @param alpha.space alpha.space is grid of values for the tuning parameter, a vector of
#' candidate values for the degree of the strength of penalty.
#' @param cv.fold cv.fold is the fold of cross validation. It should be a integer.
#' @param second.moment second.moment is a bool variable to judge whether you will add second moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.
#' @param third.moment third.moment is a bool variable to judge whether you will add third moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.
#' @param interact interact is a bool variable to judge whether you will add interaction moment conditions in balancing. It should take values in {TRUE, FALSE}. Default: TRUE.
#'
#' @return a HRB object with the following attributes:
#' \itemize{
#' \item{AT:}{ the estimate of average treatment effect in group1 (i.e, \eqn{E(Y(group1))}).}
#' \item{weight:}{ the estimated univariate balancing weight.}
#' \item{GMIM:}{ Generalized Multivariate Imbalance Measure that defines in our paper (Dai, Y., & Yan, Y. (2022).).}
#' \item{MAE:}{ Mean Absolute Error that defines in (Xu, Y., & Yang, E. (2021)).}
#' \item{alpha:}{ the tuning parameter we choose.}
#' }
#'
#' @rdname HRB
#' @references Xu, Y., & Yang, E. (2021). Hierarchically Regularized Entropy Balancing. \emph{Political Analysis}, \strong{forthcoming}, \doi{10.2139/ssrn.3807620}.
#' @references Dai, Y., & Yan, Y. (2022). Mahalanobis balancing: a multivariate perspective on approximate covariate balancing. arXiv preprint arXiv:2204.13439.
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
#' result1 <- HRB(covariate = X, treat = treat, group1 = 1, outcome = Y)
#' result2 <- HRB(covariate = X, treat = treat, group1 = 0, outcome = Y)
#'
#' # an estimate of ATE
#' result1$AT - result2$AT
#'
#' # estimating ATC
#' result3 <- HRB(covariate = X, treat = treat, group1 = 1, group2 = 0, outcome = Y)
#'
#' # an estimate of ATC
#' result3$AT - mean(data$Y[data$Tr == 0])
HRB <- function(covariate, treat, group1, group2 = NULL, outcome,
                alpha.space = c(1e-1, 0.05, 1e-2, 0.005, 1e-3),
                iterations = 1000, convergence = 1e-8, cv.fold = 5,
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
    for (k in 1:ncol(covariate)) {
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
  covariate <- scale(covariate) # Standardize covariate
  dimension <- dim(covariate)[2] # Calculate the dimension of covariate
  sample.size <- dim(covariate)[1] # Calculate the sample size of covariate
  group1.number <- sum(treat == group1) # Calculate the number of individuals in group1
  xMB <- CholMB(covariate = covariate, treat = treat, group1 = group1, group2 = group2, method = "MB") # Calculate the scaled covariate using diagonal matrix

  # Set up space for saving weight and GMIM
  weight.space <- matrix(0, nrow = group1.number, ncol = length(alpha.space))
  beta.space <- matrix(0, nrow = dimension, ncol = length(alpha.space))
  GMIM <- matrix(0, nrow = cv.fold, ncol = length(alpha.space))

  # Select the individual in group1
  x <- t(covariate[treat == group1, ])

  # Calculate the target moment condition
  if (is.null(group2)) {
    M <- colMeans(covariate)
  } else {
    M <- colMeans(covariate[treat == group2, ])
  }

  # Calculate weight for each delta in alpha.space
  fold_ids <- sample(1:group1.number) %% cv.fold + 1 # generate id for cross validation
  MAE <- matrix(NA, length(alpha.space), cv.fold) # Set up space to save mean absolute error for each cross validation
  for (i in 1:length(alpha.space)) {
    for (j in 1:cv.fold) {
      eps <- 1000 # initialize eps
      beta <- as.vector(rep(0, dimension)) # initialize beta
      x.train <- x[, fold_ids != j] # initialize training data
      x.test <- x[, fold_ids == j] # initialize test data
      Opti.func <- function(beta) {
        log(sum(exp(-crossprod(x.train, beta)))) + crossprod(M, beta) + alpha.space[i] * norm(beta, type = c("2"))^2 # loss function
      }
      # Use BFGS to optimize loss function
      beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = convergence, maxit = iterations))$par

      # Bootstrap to select tuning parameter
      weight.test <- exp(-crossprod(x.test, beta))
      weight.test <- weight.test / sum(weight.test)

      # Calculate mean absolute error for each cross validation
      MAE[j, i] <- mean(abs(tcrossprod(t(weight.test), x.test)))
    }
  }

  # calculate the AT
  meanMAE <- colMeans(MAE)
  alpha <- alpha.space[which.min(meanMAE)]

  # Using selected tuning parameter to calculate weights
  eps <- 1000 # initialize eps
  beta <- as.vector(rep(0, dimension)) # initialize beta
  Opti.func <- function(beta) {
    log(sum(exp(-crossprod(x, beta)))) + crossprod(M, beta) + alpha * norm(beta, type = c("2"))^2 # loss function
  }
  beta <- optim(par = beta, fn = Opti.func, method = "BFGS", control = list(abstol = 0, maxit = iterations))$par # Use BFGS to optimize loss function
  weight <- exp(-crossprod(x, beta)) # Calculate weight
  weight <- weight / sum(weight) # Standardize weight
  GMIM <- sum(tcrossprod(t(weight), xMB)^2) # Calculate Generalized Imbalanced Measure
  AT <- crossprod(weight, outcome[treat == group1]) # Calculate average treatment effect
  result <- list(AT = AT, weight = as.vector(weight), GMIM = GMIM, alpha = alpha, parameter = as.vector(beta))
  return(result)
}
