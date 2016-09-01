#' Print function for madr.mcmc class
#'
#' This function prints results from madr.mcmc class
#'
#' @export
#'
#' @param res madr.mcmc object
#'

print.madr.mcmc <- function(res){
   cat(paste("\nMA-DR with tau=",0))
   cat(" and two-stage procedure for calculating weights")
   cat(paste("\n\nEstimate:",res$madr))
   cat("\n\nPropensity score model inclusion probabilities:\n")
   print(round(res$weight.ps,3))
   cat("\nOutcome model inclusion probabilities:\n")
   print(round(res$weight.om,3))
}
