#' Model averaged double robust estimate with enumeration of all possible models (linear terms only)
#'
#' This function enumerates all possible models and estimates a model averaged double robust estimate
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment indicator (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param tau scalar value for the prior model dependence (1 is an independent prior)
#' @param two.stage indicator if the two-stage procedure for calculating the model weights should be used
#'
#' @export
#' @return A object of class madr.enumerate. The object contains the following named components:
#'	\item{out}{a matrix that contains the BIC and estimated treatment from each outcome model}
#'	\item{ps}{a matrix that contains the BIC from each propensity score model}
#'	\item{dr}{a matrix that contains the model-specific double robust estimates}
#'	\item{U.names}{the column names of U}
#'

madr.enumerate <- function(Y,X,U,W=NULL,tau=1,two.stage=F){
   out = OM.MA.enumerate(Y=Y,X=X,U=U,W=W)
   out. = out$out.table

   ps = PS.MA.enumerate(X,U,W=W)
   ps. = matrix(ps$out.table,ncol=2)

   dr = matrix(0,nrow(out.),nrow(ps.))
   for (i in 1:nrow(out.)){
      for (j in 1:nrow(ps.)){
         dr[i,j] = mean((Y-out[['dict']][[out.[i,1]]][['predicted']])*X/ps[['dict']][[ps.[j,1]]][['PS']]-(Y-out[['dict']][[out.[i,1]]][['predicted']])*(1-X)/(1-ps[['dict']][[ps.[j,1]]][['PS']]))  + out[['dict']][[out.[i,1]]][['beta']]
      }
   }
   colnames(dr) = paste("ps",ps.[,1],sep="")
   rownames(dr) = paste("om",out.[,1],sep="")

   res = list(out=out.,ps=ps.,dr=dr,U.names=colnames(U),tau=tau,two.stage=two.stage)
   class(res) = "madr.enumerate"
   return(res)
}
