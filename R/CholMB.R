#' @title CholMB
#' @description Function to decompose weighting matrix.
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
#' @param covariate covariate.
#' @param treat treatment indicator vector.
#' @param group1 see Details.
#' @param group2 see Details. Default: NA.
#' @param method a string that takes values in {"MB", "MB2"}, Default: 'MB'.
#' @param dimension the dimension of covariate.
#' @rdname CholMB
CholMB <- function(covariate, treat, group1, group2 = NA, method = "MB") {

  # Initiazing input
  covariate <- as.matrix(covariate) # Matrixization x
  dimension <- ncol(covariate) # Calculate the dimension of covariate
  treat <- as.vector(treat) # Vectorization treat
  covariate.group1 <- as.matrix(covariate[treat == group1, ]) # Select the covariate to be balanced
  group.number <- sum(treat == group1) # Calculate the number of individual in treatment group

  if (method == "MB") {
    meanx <- colMeans(covariate)
    var.group <- sqrt(colMeans((covariate - meanx)^2))
    group.cov <- diag(var.group)
  }
  if (method == "MB2") {
    group.cov <- cov(covariate)
  }

  mean.target <- matrix(NA, 1, dimension)
  if (is.null(group2)) {
      mean.target <- colMeans(covariate)
      covariate.group1 <- covariate.group1 - matrix(mean.target, nrow(covariate.group1), dimension, byrow = T)
  } else {
      mean.target <- colMeans(covariate[treat == group2, ])
      covariate.group1 <- covariate.group1 - matrix(mean.target, nrow(covariate.group1), dimension, byrow = T)
  }

  K <- solve(group.cov)
  Q <- chol(K)
  x <- tcrossprod(Q, covariate.group1)

  return(x)
}
