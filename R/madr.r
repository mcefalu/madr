#' Calculate model averaged double robust estimate
#'
#' This function estimates a model averaged double robust estimate.
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param M the number of MCMC iteration
#' @param cut cumulative probability of models to be retained for improved computational efficiency (1 retains all visited models)
#' @param enumerate indicator if all possible models should be enumerated (default: FALSE)
#' @param tau scalar value for the prior model dependence (1 is an independent prior; defaults to 0)
#' @param two.stage indicator if the two-stage procedure for calculating the model weights should be used (defaults to TRUE)
#'
#' @export
#' @return A list. The list contains the following named components:
#'	\item{madr}{the model averaged double robust estimate}
#'	\item{weight.ps}{a vector that contains the inclusion probability of each covariate in the propensity score model}
#'	\item{weight.om}{a vector that contains the inclusion probability of each covariate in the outcome model}
#'
madr <- function(Y,X,U,W=NULL,M=1000,cut=.95,enumerate=F,tau=NULL,two.stage=NULL){
   if (enumerate){
      if (is.null(tau)){
         tau = 0
      }
      if (is.null(two.stage)){
         two.stage = TRUE
      }
      res = madr.enumerate(Y=Y,X=X,U=U,W=W,tau=tau,two.stage=two.stage)
   }else{
      if (!is.null(tau) | !is.null(two.stage)){
         warning("\n\nUser specified values of tau and two.stage are ignored when enumerate=FALSE.\ntau is set to 0 and two.stage is set to TRUE.")
         tau = 0
         two.stage=TRUE
      }
      res <- madr.mcmc(Y=Y,X=X,U=U,W=W,M=M,cut=cut)
   }
   return(res)
}

