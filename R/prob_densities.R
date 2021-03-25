

.fAD.log_dnorm_core <- function(x,mean,sd,log){
  lp <- -0.5*square((x-mean)/sd) - 0.918938533204673 - log(sd)
  if(log){return(lp)} else { return(exp(lp))}
}


#' @export
setMethod("dnorm",c("ADtype","numeric","numeric","ANY"),
          function(x,mean=0,sd=1,log=FALSE){
            return(.fAD.log_dnorm_core(x,mean,sd,log))
          })

#' @export
setMethod("dnorm",c("ANY","ADtype","numeric","ANY"),
          function(x,mean=0,sd=1,log=FALSE){
            return(.fAD.log_dnorm_core(x,mean,sd,log))
          })

#' @export
setMethod("dnorm",c("ANY","ANY","ADtype","ANY"),
          function(x,mean=0,sd=1,log=FALSE){
            return(.fAD.log_dnorm_core(x,mean,sd,log))
          })


#.fAD.log_dgamma_core <- function(x,shape,rate=1,scale=1/rate,log=FALSE){
#  lp <- (shape-1)*log(x) - x/scale - lgamma(shape) - shape*log(scale)
#  if(log){return(lp)} else { return(exp(lp))}
#} 

# #' @export
#setMethod("dgamma",c("fAD","numeric","numeric","numeric"),
#          function(x,shape,rate=1,scale=1/rate,log=FALSE){
#            return(.fAD.log_dgamma_core(x,mean,sd,log))
#          })

# #' @export
# setMethod("dgamma",c("ANY","fAD","numeric","numeric"),
#           function(x,shape,rate=1,scale=1/rate,log=FALSE){
#             return(.fAD.log_dgamma_core(x,mean,sd,log))
#           })

# #' @export
# setMethod("dgamma",c("ANY","ANY","fAD","fAD"),
#           function(x,shape,rate=1,scale=1/rate,log=FALSE){
#             return(.fAD.log_dgamma_core(x,mean,sd,log))
#           })




