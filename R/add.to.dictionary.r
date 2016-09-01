#' Worker function that fits propensity score models
#'
#' This function fits propensity score models and saves necessary information
#'
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param alpha vector of inclusion indicators (which columns of U) to included in the propensity score model
#' @export
#' @return A list. The list contains the following named components:
#'	\item{out}{a list that contains the BIC and estimated propensity scores from propensity score models}
#'

add.to.dictionary <- function(X,U,W,alpha){
   out = list()
   if (sum(alpha)==0){
      model = glm(X~W-1,family='binomial')
   }else{
      model = glm(X~U[,alpha==1]+W-1,family='binomial')
   }
   out[['BIC']] = model$aic - 2*(sum(alpha)+1) + log(length(X))*(sum(alpha)+1)
   out[['PS']] = model$fitted
   return(out)
}
