#' Enumerates all possible outcome models (linear terms only)
#'
#' This function enumerates and fits all possible outcome models
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment indicator (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#'
#' @export
#' @return A list. The listcontains the following named components:
#'	\item{dict}{a list that contains the BIC, predicted values, and estimated treatment effect from each outcome model}
#'	\item{out.table}{a matrix that contains the BIC and estimated treatment effect from each outcome model}
#'

OM.MA.enumerate = function(Y,X,U,W=NULL){
   # add intercept to W
   W = cbind(rep(1,length(X)),W)

   # dimension of U
   k=ncol(U)
   # create all possible combinations of U
   alpha=expand.grid(as.data.frame(matrix(rep(0:1,k),ncol=k)))
   out.table = cbind(apply(alpha,1,function(x) paste(x,sep='',collapse='')),"","")
   # this is a dictionary that holds model fit info
   dict=list()
   # cycle through outcome models
   for (i in 1:nrow(out.table)){
      index = out.table[i,1]
      # fit the outcome model
      dict[[index]] = add.to.dictionary.outcome(Y=Y,X=X,U=U,W=W,alpha=alpha[i,])
      # save BIC and estimate
      out.table[i,] = c(index,dict[[index]][['BIC']],dict[[index]][['beta']])
   }
   return(list(dict=dict,out.table=out.table))
}
