#' Print function for summary.madr.enumerate class
#'
#' This function prints results from summary.madr.enumerate class
#'
#' @export
#'
#' @param res summary.madr.enumerate object
#'

print.summary.madr.enumerate <- function(res){
   cat(paste("\nMA-DR with tau=",res$tau))
   if (res$two.stage){
      cat(" and two-stage procedure for calculating weights")
   }else{
   }
   cat(paste("\n\nEstimate:",res$madr))
   cat("\n\nPropensity score model inclusion probabilities:\n")
   print(round(res$weight.ps,3))
   cat("\nOutcome model inclusion probabilities:\n")
   print(round(res$weight.om,3))
}
