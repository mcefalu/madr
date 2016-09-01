##############################
## functions for Biometrics ##
#
#  Author: Matthew Cefalu
#  Paper: A model averaged double robust estimate
#  Date last updated: 09/01/16

# expit function
expit <- function(x){
   exp(x)/(1+exp(x))
}

## enumerates all possible outcome models (linear terms only)
#
# Y = the outcome 
# X = the treatment
# U = the covariates
# W = covariates forced into the model
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

## enumerates all possible propensity score models (linear terms only)
#
# X = the treatment
# U = the covariates
# W = covariates forced into model
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

## main function, when models can be enumerated
# 
# Y = outcome
# X = treatment
# U = covariates
madr.enumerate <- function(Y,X,U,W=NULL){
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
   res = list(out=out.,ps=ps.,dr=dr,U.names=colnames(U))
   class(res) = "madr.enumerate"
   return(res)
}


## madr.enumerate post processing
summary.madr.enumerate <- function(out,tau=1,two.stage=F){
   om = out$out
   ps = out$ps
   dr = out$dr
   prior = matrix(NA,nrow(ps),nrow(om))
   rownames(prior) = paste("ps",ps[,1],sep="")
   colnames(prior) = paste("om",om[,1],sep="")
   
   for (i in 1:nrow(ps)){
      ps.mod = as.numeric(unlist(strsplit(ps[i,1],split='')))
      for (j in 1:nrow(om)){
         om.mod = as.numeric(unlist(strsplit(om[j,1],split='')))
         if (any((om.mod-ps.mod)<0)){
            prior[i,j] = tau  
         }else{
            prior[i,j] = 1
         }  
      }
   }
   # true dep prior
   a = as.numeric(om[,2])
   a = a-min(a)
   b = as.numeric(ps[,2])
   b = b-min(b)
   
   # names
   names(b) = rownames(prior)
   names(a) = colnames(prior)
   
   if (!two.stage){
      post.prob = t(t(prior)*exp(-a/2))*exp(-b/2)
      post.prob = post.prob/sum(post.prob)
   }else{
      # 2stage
      post.prob = apply(prior*exp(-b/2),2,function(x) x/sum(x))
      post.prob = t(t(post.prob)*(exp(-a/2)/sum(exp(-a/2))))
   }
   
   # ma-dr
   madr = sum(t(post.prob)*dr)
   ps = rowSums(post.prob)
   om = colSums(post.prob)
   K = length(strsplit(gsub("ps","",names(ps)[1]),split='')[[1]])
   ps.prob = numeric(K)
   om.prob = numeric(K)
   names(ps.prob) = out$U.names
   names(om.prob) = out$U.names
   
   for (i in 1:length(ps)){
      index = which(as.numeric(strsplit(gsub("ps","",names(ps)[i]),split='')[[1]])==1)
      ps.prob[index] = ps.prob[index] + ps[i]
      index = which(as.numeric(strsplit(gsub("om","",names(om)[i]),split='')[[1]])==1)
      om.prob[index] = om.prob[index] + om[i]
   }

   return(list(madr=madr,weight.ps=ps.prob,weight.om=om.prob))
}

## function to calculate model probabilities for the propensity score model
#
# outputs a dictionary that includes the BIC and the propensity score for each model visited
# 
# inputs: 
# X is treatment
# U are confounders
# W are covariates forced into the models
# alpha is a vector of inclusion indicators to start the algorithm at, can be left NULL
# master.index is the indexes of the confounders to be considered for inclusion in the PS model
#      if null, we just fit the intercept only model
#      this is necessary due to the strict prior distribution on the model space
# master.dict is a dictionary keepting information from previously visited models
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

## function that calculates model probabilities for the outcome models
#
# could be combine with PS.MA to create one general function
# Y is the outcome
# X is the treatment
# U are confounders
# W are covariates forced into the models
# M is number of MCMC
# alpha is a starting vector for inclusion indicators -- can be left null
# binary is indicator if outcome is binary
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


## adds to the PS model dictionary...fits the model and saves the BIC and propensity scores
#
# X is treatment
# U is confounders
# W are forced into the model
# alpha is vector of inclusion indicators for the confounders
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

## fits the outcome models and save useful information
#
# saves BIC, estimated beta, and predicted values...need for the DR estimator
# Y is the outcome
# X is the treatment
# U are confounders
# W are covariates forced into the models
# alpha are inclusion indicators
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

## function to convert BIC to probs
#
bic.to.prob <- function(bic){
   bf = exp(-bic/2)
   return(bf/sum(bf))
}

## main function 
#
# Y is outcome, 
# X is treatment
# U are confounders
# W are covariates forced into the models
# M is the length of MCMC
# cut is a cutoff for the cumulative posterior model probabilities.
#    cut throws out low probabilities models s.t. cumsum(probs) < cut
#    this improves computational efficiency.
madr <- function(Y,X,U,W=NULL,M=1000,cut=.95){
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
   return(res)
}

