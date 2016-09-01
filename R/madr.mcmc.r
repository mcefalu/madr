#' Calculate model averaged double robust estimate using a pseudo-MC3 algorithm
#'
#' This function uses a pseudo-MC3 algorithm to search the model space, then estimate a model averaged double robust estimate using the two-stage procedure for estimating model weights with tau=0.
#'
#' @param Y vector of the outcome
#' @param X vector of the treatment (0/1)
#' @param U matrix of covariates to be considered for inclusion/exclusion
#' @param W matrix of covariates that will be included in all models (optional)
#' @param M the number of MCMC iteration
#' @param cut cumulative probability of models to be retained for improved computational efficiency (1 retains all visited models)
#' @export
#' @return A list. The list contains the following named components:
#'	\item{madr}{the model averaged double robust estimate}
#'	\item{weight.ps}{a vector that contains the inclusion probability of each covariate in the propensity score model}
#'	\item{weight.om}{a vector that contains the inclusion probability of each covariate in the outcome model}
#'
madr.mcmc <- function(Y,X,U,W=NULL,M=1000,cut=.95){
   # first find the model probabilities for the outcome models
   out = OM.MA(Y=Y,X=X,U=U,W=W,M=M)
   # extract some information that may be useful
   out. = out$out.table
   # shift the BIC
   out.[,2] = as.numeric(out.[,2])-min(as.numeric(out.[,2]))
   # convert BIC to probs
   out. = out.[order(bic.to.prob(as.numeric(out.[,2])),decreasing=T),]
   out. = cbind(out.,bic.to.prob(as.numeric(out.[,2])))
   # throw out some models with low probs -- need to normalize the probs after the cut so that they sum to 1
   if (cut<1){
      cut.index = min(which(cumsum(out.[,4])>cut))
      out.weights = as.numeric(out.[,4])/cumsum(out.[,4])[cut.index] #cut
   }else{
      # keep all models
      cut.index = length(as.numeric(out.[,4]))
      out.weights = as.numeric(out.[,4])/cumsum(out.[,4])[cut.index]
   }
   # create a dictionary for the propensity score models
   # this avoids refitting a model we already visited
   master.dict = list()
   # keep track of the ma.dr as we go
   ma.dr = 0
   # save PS model weights
   ps.out = numeric(ncol(U))
   om.out = numeric(ncol(U))
   # loop over all outcome models that were visited -- the prior restricts the PS model space based on the outcome model
   # under independent prior, this loop is not needed. just fit PS models independently.
   for (i in 1:cut.index){
      if ((i%%5)==0){
         print(paste('Outcome model',i,'of',cut.index))
      }
      index = NULL
      # check that the outcome model includes at least one confounder
      if (any(as.numeric(strsplit(out.[i,1],split='')[[1]])==1)){
         # if it does, grab the indexes for the confounders included in the outcome model
         index = which(as.numeric(strsplit(out.[i,1],split='')[[1]])==1)
      }
      # find propensity score model probabilities for this outcome model
      # one additional alteration to make the code faster is to make M here a function of the length of index
      ps = PS.MA(X,U,W=W,M=M,master.index=index,master.dict=master.dict)
      # save the dictionary for next iteration
      master.dict = ps$dict
      # save other things that may be useful
      ps. = matrix(ps$out.table,ncol=2)
      # rescale BIC
      ps.[,2] = as.numeric(ps.[,2])-min(as.numeric(ps.[,2]))
      # convert BIC to probs
      ps. = matrix(ps.[order(bic.to.prob(as.numeric(ps.[,2])),decreasing=T),],ncol=2)
      ps. = cbind(ps.,bic.to.prob(as.numeric(ps.[,2])))
      # do the cut of the model probs
      # this step only speeds up computation slightly
      if (cut<1){
         ps.index = min(which(cumsum(ps.[,3])>cut))
         ps.weights = as.numeric(ps.[,3])/cumsum(ps.[,3])[ps.index]#/cut
      }else{
         ps.index = length(ps.[,3])
         ps.weights = as.numeric(ps.[,3])/cumsum(ps.[,3])[ps.index]#/cut
      }
      # loop over the ps models that were visited
      for (j in 1:ps.index){
         # calculated DR estimator for outcome model indexed by i and PS model indexed by j
         dr = mean((Y-out[['dict']][[out.[i,1]]][['predicted']])*X/ps[['dict']][[ps.[j,1]]][['PS']]-(Y-out[['dict']][[out.[i,1]]][['predicted']])*(1-X)/(1-ps[['dict']][[ps.[j,1]]][['PS']]))  + out[['dict']][[out.[i,1]]][['beta']]
         # now weight it based on its posterior model probability
         ma.dr = ma.dr + out.weights[i]*ps.weights[j]*dr
         # save weights by variable
         ps.out[which(as.numeric(strsplit(ps.[j,1],split='')[[1]])==1)] = ps.weights[j]*out.weights[i] + ps.out[which(as.numeric(strsplit(ps.[j,1],split='')[[1]])==1)]
      }
      # save outcome model weights
      om.out[which(as.numeric(strsplit(out.[i,1],split='')[[1]])==1)] = out.weights[i] + om.out[which(as.numeric(strsplit(out.[i,1],split='')[[1]])==1)]
   }
   names(ps.out) = colnames(U)
   names(om.out) = colnames(U)
   res = list(madr=ma.dr,weight.ps=ps.out,weight.om=om.out)
   class(res) = "madr.mcmc"
   return(res)
}

