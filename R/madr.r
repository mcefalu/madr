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
#' @importFrom stats coef glm rbinom runif
#'
#' @examples
#' set.seed(122)
#' ## generate data
#' n = 100 # number of observations
#' k = 4   # number of covariates
#' U = matrix(rnorm(n*k),n,k)
#' colnames(U) = paste0("U",1:k)
#' A = rbinom(n,1,expit(-1+.5*rowSums(U)))
#' Y = rnorm(n,1+A+.25*rowSums(U))
#'
#' ## A is confounded -- true effect is 1
#' lm(Y~A)
#'
#' ## fit ma-dr -- can enumerate models if k isnt too big
#' res = madr(Y=Y,X=A,U=U,enumerate=TRUE,tau=1,two.stage=FALSE) # independent prior
#' res
#'
#' res = madr(Y=Y,X=A,U=U,enumerate=TRUE,tau=0,two.stage=TRUE) # tau=0 and using two-stage weights
#' res
#'
#' ## no need to refit madr each time when enumerating -- use summarize and specify different taus
#' summary(res,tau=1,two.stage=FALSE) # independent prior
#' summary(res,tau=0,two.stage=FALSE)
#' summary(res,tau=0,two.stage=TRUE) # two-stage procedure for calculating weights
#'
#' ## use mcmc instead of enumerating (the default)
#' madr(Y=Y,X=A,U=U,M=1000,cut=1) #should approximate tau=0 and two.stage=TRUE

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

