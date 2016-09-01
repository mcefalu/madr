#' Calculate model probabilities for the propensity score model using a pseudo-MC3 algorithm
#'
#' This function uses a pseudo-MC3 algorithm to search the propensity score model space.
#'
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param M the number of MCMC iteration
#' @param alpha vector of inclusion indicators (which columns of U) to start MCMC algorithm (optional)
#' @param master.index indexes which columns of U should be considered for inclusion in the propensity score model (optional)
#' @param master.dict list containing information from previous propensity score model fits (optional)
#' @export
#' @return A list. The list contains the following named components:
#'	\item{dict}{a list that contains the BIC and estimated propensity scores from propensity score models}
#'	\item{alpha}{the last model visited by the algorithm}
#'	\item{out.table}{a matrix that contains the BIC from each propensity score model}
#'

PS.MA = function(X,U,W=NULL,M=1000,alpha=NULL,master.index=NULL,master.dict=list()){
   # add intercept to W
   W = cbind(rep(1,length(X)),W)

   # uses a psuedo mc3 algorithm to jump model
   k=ncol(U)

   # if the master.index is null (which would mean there are no confounders included in the outcome model), we just fit the intercept only model
   if (is.null(master.index)){
      # the index of the dictionary is just a string of 0/1 corresponding to whether each variable is included in the model
      # here is only 0's
      index = numeric(k)
      index = paste(index,sep='',collapse='')
      dict = master.dict
      dict[[index]]=list()
      # fit / add the null model to the dictionary
      model = glm(X~W-1,family='binomial')
      dict[[index]][['BIC']] = model$aic - 2*(1+1) + log(length(X))*(1+1)
      dict[[index]][['PS']] = model$fitted
      out.table = c(index,dict[[index]][['BIC']])
      return(list(dict=dict,out.table=out.table))
   }

   # this should never be visited...
   if (is.null(master.index)){
      master.index = 1:ncol(U)
   }

   # restrict U to only include those confounder under consideration -- i.e. those in master.index -- i.e. those in the outcome model
   U = as.matrix(U[,master.index])

   # check if we know where we want to start the set of inclusion indicators for the confounders
   if (is.null(alpha)){
      # generate random starting vector of inclusion indicators
      alpha = rbinom(ncol(U),1,.5)
      while (sum(alpha)>(nrow(U)/1.5)){
         alpha = rbinom(ncol(U),1,nrow(U)/ncol(U)/2)
      }
   }
   # list holding the model fit information and the count the model has been visited
   dict = master.dict # first index is the model
   rm(master.dict)
   # note that there is some bookkeeping here -- the indexes within the function are different than the global indexes
   # therefore, we just always make sure we use the global indexes
   # the first column of U here is not necessary the first column of U outside of the function
   index = numeric(k)
   index[alpha*master.index]=1
   index = paste(index,sep='',collapse='')
   # add.to.dictionary outputs the information we store in the dictionary
   dict[[index]] = add.to.dictionary(X=X,U=U,W=W,alpha=alpha)
   # the count is not important
   dict[[index]][['count']] =  1
   out.table = c(index,dict[[index]][['BIC']])

   for (i in 2:M){
      # first choose variable to add or remove
      index.prop = sample(length(alpha),1)
      if((sum(alpha)+1)>=nrow(U)){
         index.prop = sample(which(alpha==1),1)
      }
      # these if/else can be collapsed...but it is useful to keep separate
      if (alpha[index.prop]==0){ # add a variable!
         alpha.prop = alpha
         alpha.prop[index.prop] = 1
      }else{ #remove a variable!
         alpha.prop = alpha
         alpha.prop[index.prop] = 0
      }
      # switch an included/excluded -- but only with small probability
      if ( (any(alpha==1)) & (any(alpha==0)) ){
         if (runif(1)<.25){
            index.prop = sample(which(alpha==1),1)
            index.prop2 = sample(which(alpha==0),1)
            # note that we overwrite the previous proposed move
            alpha.prop = alpha
            alpha.prop[index.prop] = 0
            alpha.prop[index.prop2] = 1
         }
      }
      # find the index of the proposed move -- note the bookkeeping again
      index.prop = numeric(k)
      index.prop[alpha.prop*master.index]=1
      index.prop = paste(index.prop,sep='',collapse='')
      # only fit the model if we have not visited it yet! so check that it is not in the dictionary
      if (is.null(dict[[index.prop]])){
         dict[[index.prop]] = add.to.dictionary(X=X,U=U,W=W,alpha=alpha.prop)
         out.table = rbind(out.table,c(index.prop,dict[[index.prop]][['BIC']]))
      }
      # convert BIC to probs
      num = dict[[index.prop]][['BIC']]
      den = dict[[index]][['BIC']]
      den = den - num
      num = 0

      # accept/reject the move!
      if ((exp(-num/2)/(exp(-num/2)+exp(-den/2)))>runif(1)){
         alpha = alpha.prop
         index = index.prop
      }else{ }
      # count how many time we end up in this model -- again, not important
      dict[[index]][['count']] = dict[[index]][['count']] + 1
   }
   return(list(dict=dict,alpha=alpha,out.table=out.table))
}
