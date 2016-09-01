#' Enumerates all possible propensity score models (linear terms only)
#'
#' This function enumerates and fits all possible propensity score models
#'
#' @param X vector of the treatment indicator (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#'
#' @export
#' @return A list. The list contains the following named components:
#'	\item{dict}{a list that contains the BIC and estimated propensity scores from propensity score models}
#'	\item{out.table}{a matrix that contains the BIC from each propensity score model}
#'

PS.MA.enumerate = function(X,U,W=NULL){
   # add intercept to W
   W = cbind(rep(1,length(X)),W)

   # dimension of U
   k=ncol(U)
   # create all possible combinations of U
   alpha=expand.grid(as.data.frame(matrix(rep(0:1,k),ncol=k)))
   out.table = cbind(apply(alpha,1,function(x) paste(x,sep='',collapse='')),"")
   # empty dictionary to hold model fits
   dict=list()
   for (i in 1:nrow(out.table)){
      index = out.table[i,1]
      dict[[index]] = add.to.dictionary(X=X,U=U,W=W,alpha=alpha[i,])
      # save BIC
      out.table[i,] = c(index,dict[[index]][['BIC']])
   }
   return(list(dict=dict,out.table=out.table))
}
