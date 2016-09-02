#' Print function for madr.mcmc class
#'
#' This function prints results from madr.mcmc class
#'
#' @export
#'
#' @param x madr.mcmc object
#' @param ... ignored
#'

print.madr.mcmc <- function(x,...){
   cat(paste("\nMA-DR with tau=",0))
   cat(" and two-stage procedure for calculating weights")
   cat(paste("\n\nEstimate:",x$madr))
   cat("\n\nPropensity score model inclusion probabilities:\n")
   print(round(x$weight.ps,3))
   cat("\nOutcome model inclusion probabilities:\n")
   print(round(x$weight.om,3))
}
