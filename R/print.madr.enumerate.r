#' Print function for madr.enumerate class
#'
#' This function prints results from madr.enumerate class
#'
#' @export
#'
#' @param out madr.enumerate object
#'

print.madr.enumerate <- function(out){
   res = summary(out)
   cat(paste("\nMA-DR with tau=",out$tau))
   if (out$two.stage){
      cat(" and two-stage procedure for calculating weights")
   }else{
   }
   cat(paste("\n\nEstimate:",res$madr))
   cat("\n\nPropensity score model inclusion probabilities:\n")
   print(round(res$weight.ps,3))
   cat("\nOutcome model inclusion probabilities:\n")
   print(round(res$weight.om,3))
}
