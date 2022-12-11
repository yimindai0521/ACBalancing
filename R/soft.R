#' soft function
#'
#' @description proximal operator for \eqn{l_1} norm and \eqn{l_2} norm
#'
#' @param vt input for proximal operator, it should be a vector.
#' @param lambda threshold for proximal operator, it should be positive.
#' @param norm norm type, it should take values in c("square-l2", "l1").
#'
#' @return soft_value output for proximal operator
soft <- function(vt, lambda, norm = "square-l2") {
  if (norm == "square-l2") {
    normvalue <- sqrt(sum(vt^2))
    if (normvalue > lambda) {
      return((1 - lambda / normvalue) * vt) # output = (1 - lambda / ||beta||_2) * beta if ||beta||_2 >= lambda
    } else {
      return(rep(0, length(vt))) # output = 0 if lambda > ||beta||_2
    }
  }
  if (norm == "l1") {
    soft_value <- sign(vt) * abs(abs(vt) - lambda) # update = sign(a) * ||a| - lambda|
    soft_value[abs(vt) < lambda] <- 0 # update = |a| - lambda > 0
    return(soft_value)
  }
}
