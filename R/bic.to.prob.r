#' Convert BIC to model probabilities
#'
#' This function transforms BIC to model probabilities
#'
#' @param bic vector of BICs
#'
#' @export
#' @return A vector of model probabilities of the same dimension of \code{bic}
#'
bic.to.prob <- function(bic){
   bf = exp(-bic/2)
   return(bf/sum(bf))
}
