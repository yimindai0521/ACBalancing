#' soft
#'
#' @param vt
#' @param lambda
#' @param norm
#'
#' @return
#'
#' @examples
soft <- function(vt, lambda, norm = "square-l2") {
  if(norm == "square-l2"){
    normvalue <- sqrt(sum(vt^2))
    if (normvalue > lambda) {
      return((1 - lambda / normvalue) * vt)
    } else {
    return(rep(0, length(vt)))
    }
  }
  if(norm == "l1"){
    soft_value <- sign(vt) * abs(abs(vt) - lambda) # update = sign(a) * ||a| - lambda|
    soft_value[abs(vt) < lambda] <- 0 # update = |a| - lambda > 0
    return(soft_value)
  }

}
