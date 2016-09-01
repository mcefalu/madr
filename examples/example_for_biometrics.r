## Example of model averaged double robust estimation

rm(list=ls())

library(madr)

###### example #########
   set.seed(122)

   n = 100 # number of observations
   k = 4   # number of covariates

## generate data
   U = matrix(rnorm(n*k),n,k)
   colnames(U) = paste0("U",1:k)
   A = rbinom(n,1,expit(-1+.5*rowSums(U)))
   Y = rnorm(n,1+A+.25*rowSums(U))

## A is confounded -- true effect is 1
   lm(Y~A)

## fit ma-dr -- can enumerate models if k isnt too big
   res = madr(Y=Y,X=A,U=U,enumerate=T,tau=1,two.stage=F) # independent prior
   res

   res = madr(Y=Y,X=A,U=U,enumerate=T,tau=0,two.stage=T) # dependent prior with tau=0 and using two-stage weights
   res

## no need to refit madr each time when enumerating -- use summarize and specify different taus
   summary(res,tau=1,two.stage=F) # independent prior
   summary(res,tau=0,two.stage=F)
   summary(res,tau=0,two.stage=T) # two-stage procedure for calculating weights

## use mcmc instead of enumerating (the default)
## tau and two.stage are not used when enumerate=F (the default)
   madr(Y=Y,X=A,U=U,M=1000,cut=1,tau=4) # should approximate tau=0 and two.stage=T

## for fun, plot madr as a function of tau
   Tau = seq(0,1,.01)
   madr = NULL
   for (tau in Tau){
      madr = c(madr,summary(res,tau=tau)$madr)
   }
   plot(x=Tau,y=madr,type='l')

# also, histogram of all model specific DR
   hist(res$dr,main="Model-specific double robust estimates",xlab="DR")



#####################
## larger set of U ##
   set.seed(542)
   n = 100 # number of observations
   k = 50   # number of covariates
   U = matrix(rnorm(n*k),n,k)
   colnames(U) = paste0("U",1:k)

## 1-5 are related to treatment
   A = rbinom(n,1,expit(-1+rowSums(U[,1:5])))

# only 4-5 are related to outcome
   Y = rnorm(n,1+A+rowSums(U[,4:5]))

# A is confounded
   lm(Y~A)

## fit ma-dr
#
# I set M = phi*(the # of confounders), where phi is the number of times
# you want each variable selected for inclusion/exclusion (on average)
# in my simulations, i set phi to be 50
#
# A value of cut closer to 1 keeps more models -- thus takes longer to run
#
# as noted in the madr function, you can adaptively change M when fitting the PS model
# based on the # of confounders included in the outcome
   madr(Y=Y,X=A,U=U,M=20*ncol(U),cut=.9)




