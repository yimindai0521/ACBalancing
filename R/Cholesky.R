#' @title Cholesky
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
#' @param x covariates.
#' @param treat treatment indicator vector.
#' @param group1 see Details.
#' @param group2 see Details. Default: NA.
#' @param method a string that takes values in {"MB", "MB2"}, Default: 'MB'.
#' @param dimension the dimension of covariate.
#' @rdname Cholesky
#' @export
Cholesky <- function(x, treat, group1, group2 = NA, method = "MB", dimension) {
  data.matrix <- cbind(x, treat)
  group <- as.matrix(data.matrix[treat == group1, 1:dimension])
  group.num <- sum(treat == group1)
  group.col <- ncol(group)
  group.row <- nrow(group)
  if (is.null(group.col)) {
    group.col <- 1
    group.row <- sum(abs(group) >= 0)
  }
  group.cov <- matrix(0, group.col, group.col)

  if (method == "MB") {
    for (i in 1:group.col) {
      meanx <- sum(group[, i]) / group.row
      group.cov[i, i] <- sum((group[, i] - meanx) * (group[, i] - meanx)) / (group.num - 1)
    }
  }
  if (method == "MB2") {
    for (i in 1:group.col) {
      for (j in 1:group.col) {
        meanx <- sum(group[, i]) / group.row
        meany <- sum(group[, j]) / group.row
        group.cov[i, j] <- sum((group[, i] - meanx) * (group[, j] - meany)) / (group.num - 1)
      }
    }
  }

  if (is.na(group2)) {
    mean.pop <- matrix(NA, 1, dimension)
    if (dimension == 1) {
      mean.pop <- mean(x)
      group.se <- data.matrix[(treat == group1), 1] - mean.pop
    } else {
      mean.pop <- apply(data.matrix[, 1:dimension], 2, mean)
      group.se <- data.matrix[(treat == group1), 1:dimension]
      for (i in 1:sum(treat == group1)) {
        group.se[i, ] <- group.se[i, ] - mean.pop
      }
    }
  } else {
    mean.pop <- matrix(NA, 1, dimension)
    if (dimension == 1) {
      mean.pop <- mean(x[treat == group2])
      group.se <- data.matrix[(treat == group1), 1] - mean.pop
    } else {
      mean.pop <- apply(data.matrix[treat == group2, 1:dimension], 2, mean)
      group.se <- data.matrix[(treat == group1), 1:dimension]
      for (i in 1:sum(treat == group1)) {
        group.se[i, ] <- group.se[i, ] - mean.pop
      }
    }
  }

  K <- solve(group.cov)
  Q <- chol(K)
  x <- Q %*% t(group.se)

  return(x)
}
