#' Expit (inverse logit) function
#'
#' This function transforms the input using the expit function
#'
#' @param x vector of values to apply the expit function
#'
#' @export
#' @return A vector of the same dimension of \code{x}
#'
#'
expit <- function(x){
   exp(x)/(1+exp(x))
}
