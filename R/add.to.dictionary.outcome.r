#' Worker function that fits outcome models
#'
#' This function fits outcome models and saves necessary information
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param alpha vector of inclusion indicators (which columns of U) to included in the propensity score model
#' @param binary indicates if the outcome is binary
#' @export
#' @return A list. The list contains the following named components:
#'	\item{out}{a list that contains the BIC, predicted values, and estimated treatment effect from each outcome model}
#'
add.to.dictionary.outcome <- function(Y,X,U,W,alpha,binary=F){
   if (!binary){
      out = list()
      if (sum(alpha)==0){
         model = glm(Y~X+W-1,family='gaussian')
      }else{
         model = glm(Y~X+U[,alpha==1]+W-1,family='gaussian')
      }
      out[['BIC']] = model$aic - 2*(sum(alpha)+1) + log(length(X))*(sum(alpha)+1)
      out[['beta']] = coef(model)[1]
      out[['predicted']] = model$fitted
      return(out)
   }else{
      out = list()
      if (sum(alpha)==0){
         model = glm(Y~X+W-1,family='binomial')
      }else{
         model = glm(Y~X+U[,alpha==1]+W-1,family='binomial')
      }
      out[['BIC']] = model$aic - 2*(sum(alpha)+1) + log(length(X))*(sum(alpha)+1)
      out[['beta']] = mean(expit(cbind(1,U[,alpha==1],W)%*%coef(model))-expit(cbind(0,U[,alpha==1],W)%*%coef(model)))
      out[['predicted']] = model$fitted
      return(out)
   }
}
