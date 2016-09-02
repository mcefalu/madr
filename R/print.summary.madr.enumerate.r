#' Print function for summary.madr.enumerate class
#'
#' This function prints results from summary.madr.enumerate class
#'
#' @export
#'
#' @param x summary.madr.enumerate object
#' @param ... ignored
#'

print.summary.madr.enumerate <- function(x,...){
   cat(paste("\nMA-DR with tau=",x$tau))
   if (x$two.stage){
      cat(" and two-stage procedure for calculating weights")
   }else{
   }
   cat(paste("\n\nEstimate:",x$madr))
   cat("\n\nPropensity score model inclusion probabilities:\n")
   print(round(x$weight.ps,3))
   cat("\nOutcome model inclusion probabilities:\n")
   print(round(x$weight.om,3))
}
