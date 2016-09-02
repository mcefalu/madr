#' Provides model averaged double robust estimate for different values of tau
#'
#' This function estimates model averaged double robust estimate for different values of tau using a madr.enumerate object
#'
#' @param object madr.enumerate object
#' @param tau scalar value for the prior model dependence (1 is an independent prior; defaults to value used in madr.enumerate)
#' @param two.stage indicator if the two-stage procedure for calculating the model weights should be used (defaults to value used in madr.enumerate)
#' @param ... ignored
#'
#' @export
#' @return A list. The list contains the following named components:
#'	\item{madr}{the model averaged double robust estimate}
#'	\item{weight.ps}{a vector that contains the inclusion probability of each covariate in the propensity score model}
#'	\item{weight.om}{a vector that contains the inclusion probability of each covariate in the outcome model}
#'	\item{tau}{value of tau used in estimation}
#'	\item{two.stage}{indicator if the two-stage procedure for calculating the model weights was used}
#'

summary.madr.enumerate <- function(object,tau=NULL,two.stage=NULL,...){
   if (is.null(tau)){
      tau = object$tau
   }
   if (is.null(two.stage)){
      two.stage = object$two.stage
   }
   om = object$out
   ps = object$ps
   dr = object$dr
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
   names(ps.prob) = object$U.names
   names(om.prob) = object$U.names

   for (i in 1:length(ps)){
      index = which(as.numeric(strsplit(gsub("ps","",names(ps)[i]),split='')[[1]])==1)
      ps.prob[index] = ps.prob[index] + ps[i]
      index = which(as.numeric(strsplit(gsub("om","",names(om)[i]),split='')[[1]])==1)
      om.prob[index] = om.prob[index] + om[i]
   }
   res = list(madr=madr,weight.ps=ps.prob,weight.om=om.prob,tau=tau,two.stage=two.stage)
   class(res) = "summary.madr.enumerate"
   return(res)
}
