#' Calculate model probabilities for the outcome models using a pseudo-MC3 algorithm
#'
#' This function uses a pseudo-MC3 algorithm to search the outcome model space.
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param M the number of MCMC iteration
#' @param alpha vector of inclusion indicators (which columns of U) to start MCMC algorithm (optional)
#' @param binary indicator if the outcome is binary (optional)
#' @export
#' @return A list. The list contains the following named components:
#'	\item{dict}{a list that contains the BIC, predicted values, and estimated treatment effect from each outcome model}
#'	\item{alpha}{the last model visited by the algorithm}
#'	\item{out.table}{a matrix that contains the BIC and estimated treatment effect from each outcome model}
#'
OM.MA = function(Y,X,U,W=NULL,M=1000,alpha=NULL,binary=F){
   # add intercept to W
   W = cbind(rep(1,length(X)),W)

   # uses an pseudo-mc3 algorithm to jump models
   # create starting value for inclusion indicators if one is not provided
   if (is.null(alpha)){
      alpha = rbinom(ncol(U),1,.5)
      while (sum(alpha)>(nrow(U)/1.5)){
         alpha = rbinom(ncol(U),1,nrow(U)/ncol(U)/2)
      }
   }
   # list holding the model fit information and the count the model has been visited
   dict = list() # first index is the model
   index = paste(alpha,sep='',collapse='')
   dict[[index]] = add.to.dictionary.outcome(Y=Y,X=X,U=U,W=W,alpha=alpha,binary=binary)
   dict[[index]][['count']] =  1
   out.table = c(index,dict[[index]][['BIC']],dict[[index]][['beta']])

   for (i in 2:M){
      # first choose variable to add or remove
      index.prop = sample(length(alpha),1)
      if((sum(alpha)+1)>=nrow(U)){
         index.prop = sample(which(alpha==1),1)
      }
      if (alpha[index.prop]==0){ # add a variable!
         alpha.prop = alpha
         alpha.prop[index.prop] = 1
      }else{ # remove a variable
         alpha.prop = alpha
         alpha.prop[index.prop] = 0
      }
      # switch a included/excluded
      if ( (any(alpha==1)) & (any(alpha==0)) ){
         if (runif(1)<.25){
            index.prop = sample(which(alpha==1),1)
            index.prop2 = sample(which(alpha==0),1)
            # overwrite previous proposed alpha
            alpha.prop = alpha
            alpha.prop[index.prop] = 0
            alpha.prop[index.prop2] = 1
         }
      }

      index.prop = paste(alpha.prop,sep='',collapse='')
      # only fit the model if it is not in the dictionary
      if (is.null(dict[[index.prop]])){
         dict[[index.prop]] = add.to.dictionary.outcome(Y=Y,X=X,U=U,W=W,alpha=alpha.prop,binary=binary)
         out.table = rbind(out.table,c(index.prop,dict[[index.prop]][['BIC']],dict[[index]][['beta']]))
      }
      # convert BIC to probs
      num = dict[[index.prop]][['BIC']]
      den = dict[[index]][['BIC']]
      den = den - num
      num = 0

      # accept or reject the move!
      if ((exp(-num/2)/(exp(-num/2)+exp(-den/2)))>runif(1)){
         # keep proposed value
         alpha = alpha.prop
         index = index.prop
      }else{ }
      # count how many time we end up in this model -- not important
      dict[[index]][['count']] = dict[[index]][['count']] + 1
   }
   return(list(dict=dict,alpha=alpha,out.table=out.table))
}
