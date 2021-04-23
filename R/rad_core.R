


.as.fAD <- function(x){
  return(new("fAD",val=x,jac=matrix(0.0,length(x),.GlobalEnv$fAD.nvars__)))
}

setMethod("show",c("fAD"),function(object){
  print(paste0("fAD object, length=",length(object@val), ", values : "),quote = FALSE)
  print(object@val)
  print("Jacobian : ",quote = FALSE)
  print(object@jac)
})

setMethod("show",c("fAD2"),function(object){
  print(paste0("fAD2 object, length=",length(object@val), ", values : "),quote = FALSE)
  print(object@val)
  print("Jacobian : ",quote = FALSE)
  print(value(object@jac))
  print("hessian : ",quote = FALSE)
  for(i in 1:length(object@val)){
    print(paste0("dim # : ",i),quote=FALSE)
    print(object@hessian[i,,])
  }
})


setMethod("as.vector",c("fAD"),function(x){
  return(x)
})
setMethod("as.vector",c("fAD2"),function(x){
  return(x)
})

setMethod("length",c("fAD"),function(x){return(length(x@val))})
setMethod("length",c("fAD2"),function(x){return(length(x@val))})


setMethod("[",c("fAD","ANY","ANY","ANY"),
          function (x, i, j, drop = TRUE){
            if(! missing(j)) stop("too many arguments to []")
            
            if(length(i)==1){
              n.jac <- matrix(x@jac[i,],1,.GlobalEnv$fAD.nvars__)
            } else {
              n.jac <- x@jac[i,]
            }
            return(new("fAD",val=x@val[i],jac=n.jac))
          })
setMethod("[",c("fAD2","ANY","ANY","ANY"),
          function (x, i, j, drop = TRUE){
            if(! missing(j)) stop("too many arguments to []")
            
            if(length(i)==1){
              n.jac <- matrix(x@jac[i,],1,.GlobalEnv$fAD.nvars__)
              n.hess <- array(x@hessian[i,,],c(1,.GlobalEnv$fAD.nvars__,.GlobalEnv$fAD.nvars__))
            } else {
              n.jac <- x@jac[i,]
              n.hess <- x@hessian[i,,]
            }
            return(new("fAD2",val=x@val[i],jac=n.jac,hessian=n.hess))
          })
# faster version of []
#' @export
setGeneric("coeff",function(x,i) standardGeneric("coeff"))
#' @export
setMethod("coeff",c("fAD","ANY"),function(x,i){
  return(new("fAD",val=x@val[i],jac=matrix(x@jac[i,],length(i))))
})

#' @export
setMethod("coeff",c("fAD2","ANY"),function(x,i){
  return(new("fAD2",val=x@val[i],jac=matrix(x@jac[i,],length(i)),
             hessian=array(x@hessian[i,,],c(length(i),.GlobalEnv$fAD.nvars__,.GlobalEnv$fAD.nvars__))))
})

#' @export
setMethod("coeff",c("numeric","ANY"),function(x,i){
  return(return(x[i]))
})

setMethod("[<-",
          signature(x = "fAD","ANY","ANY","ANY"),
          function (x, i, j, ..., value) 
          {
            
            if(class(value)=="numeric"){
              x@val[i] <- value
              x@jac[i,] <- 0.0
            } else if(class(value)=="fAD"){
              x@val[i] <- value@val
              x@jac[i,] <- value@jac
            } else {
              stop("bad right hand side in [<-")
            }
            return(x)
          })

setMethod("[<-",
          signature(x = "fAD2","ANY","ANY","ANY"),
          function (x, i, j, ..., value) 
          {
            
            if(class(value)=="numeric"){
              x@val[i] <- value
              x@jac[i,] <- 0.0
              x@hessian[i,,] <- 0.0
            } else if(class(value)=="fAD2"){
              x@val[i] <- value@val
              x@jac[i,] <- value@jac
              x@hessian[i,,] <- value@hessian
            } else {
              stop("bad right hand side in [<-")
            }
            return(x)
          })

#' @export
setGeneric("coeffRef<-",function(x,i,value) standardGeneric("coeffRef<-") )

#' @export
setMethod("coeffRef<-",c("fAD","ANY","fAD"),
          function(x,i,value){
            x@val[i] <- value@val
            x@jac[i,] <- value@jac
            return(x)
          })

#' @export
setMethod("coeffRef<-",c("fAD2","ANY","fAD2"),
          function(x,i,value){
            x@val[i] <- value@val
            x@jac[i,] <- value@jac
            x@hessian[i,,] <- value@hessian
            return(x)
          })


#
# rep-like methods
#


#' @export
setMethod("rep.int",c("fAD","ANY"),function(x,times){ 
  inds <- rep.int(1:length(x),times)
  return(new("fAD",val=x@val[inds],jac=x@jac[inds,]))
})

#' @export
setMethod("rep.int",c("fAD2","ANY"),function(x,times){ 
  inds <- rep.int(1:length(x),times)
  return(new("fAD2",val=x@val[inds],jac=x@jac[inds,],hessian=x@hessian[inds,,]))
})


#' @export
setMethod("rep_len",c("fAD","ANY"),function(x,length.out){
  inds <- rep_len(1:length(x),length.out)
  return(new("fAD",val=x@val[inds],jac=x@jac[inds,]))
})

#' @export
setMethod("rep_len",c("fAD2","ANY"),function(x,length.out){
  inds <- rep_len(1:length(x),length.out)
  return(new("fAD2",val=x@val[inds],jac=x@jac[inds,],hessian=x@hessian[inds,,]))
})




.fAD.zeros <- function(n){
  return(new("fAD",val=rep.int(0.0,n),jac=matrix(0.0,n,.GlobalEnv$fAD.nvars__)))
}

.fAD2.zeros <- function(n){
  return(new("fAD2",val=rep.int(0.0,n),jac=matrix(0.0,n,.GlobalEnv$fAD.nvars__),
             hessian=array(0.0,c(n,.GlobalEnv$fAD.nvars__,.GlobalEnv$fAD.nvars__))))
}

.ADtype.zeros <- function(n,type,type2=NULL){
  if(!is.null(type2)){
    if(any(type2==.AD_types)) type<-type2
  }
  
  if(type=="fAD"){
    return(.fAD.zeros(n))
  } else if(type=="fAD2"){
    return(.fAD2.zeros(n))
  } else {
    
    stop("bad type")
  }
}

#
# binary non-numerical operations
#


# ==
setMethod("==",c("fAD","fAD"),function(e1,e2){
  return(e1@val==e2@val)
})
setMethod("==",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val==e2@val)
})
setMethod("==",c("fAD","numeric"),function(e1,e2){
  return(e1@val==e2)
})
setMethod("==",c("fAD2","numeric"),function(e1,e2){
  return(e1@val==e2)
})
setMethod("==",c("numeric","fAD"),function(e1,e2){
  return(e1==e2@val)
})
setMethod("==",c("numeric","fAD2"),function(e1,e2){
  return(e1==e2@val)
})

# !=
setMethod("!=",c("fAD","fAD"),function(e1,e2){
  return(e1@val!=e2@val)
})
setMethod("!=",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val!=e2@val)
})
setMethod("!=",c("fAD","numeric"),function(e1,e2){
  return(e1@val!=e2)
})
setMethod("!=",c("fAD2","numeric"),function(e1,e2){
  return(e1@val!=e2)
})
setMethod("!=",c("numeric","fAD"),function(e1,e2){
  return(e1!=e2@val)
})
setMethod("!=",c("numeric","fAD2"),function(e1,e2){
  return(e1!=e2@val)
})





#>
setMethod(">",c("fAD","fAD"),function(e1,e2){
  return(e1@val>e2@val)
})
setMethod(">",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val>e2@val)
})
setMethod(">",c("fAD","numeric"),function(e1,e2){
  return(e1@val>e2)
})
setMethod(">",c("fAD2","numeric"),function(e1,e2){
  return(e1@val>e2)
})
setMethod(">",c("numeric","fAD"),function(e1,e2){
  return(e1>e2@val)
})
setMethod(">",c("numeric","fAD2"),function(e1,e2){
  return(e1>e2@val)
})


#>=
setMethod(">=",c("fAD","fAD"),function(e1,e2){
  return(e1@val>=e2@val)
})
setMethod(">=",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val>=e2@val)
})
setMethod(">=",c("fAD","numeric"),function(e1,e2){
  return(e1@val>=e2)
})
setMethod(">=",c("fAD2","numeric"),function(e1,e2){
  return(e1@val>=e2)
})
setMethod(">=",c("numeric","fAD"),function(e1,e2){
  return(e1>=e2@val)
})
setMethod(">=",c("numeric","fAD2"),function(e1,e2){
  return(e1>=e2@val)
})






# <
setMethod("<",c("fAD","fAD"),function(e1,e2){
  return(e1@val<e2@val)
})
setMethod("<",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val<e2@val)
})
setMethod("<",c("fAD","numeric"),function(e1,e2){
  return(e1@val<e2)
})
setMethod("<",c("fAD2","numeric"),function(e1,e2){
  return(e1@val<e2)
})
setMethod("<",c("numeric","fAD"),function(e1,e2){
  return(e1<e2@val)
})
setMethod("<",c("numeric","fAD2"),function(e1,e2){
  return(e1<e2@val)
})


# <=
setMethod("<=",c("fAD","fAD"),function(e1,e2){
  return(e1@val<=e2@val)
})
setMethod("<=",c("fAD2","fAD2"),function(e1,e2){
  return(e1@val<=e2@val)
})
setMethod("<=",c("fAD","numeric"),function(e1,e2){
  return(e1@val<=e2)
})
setMethod("<=",c("fAD2","numeric"),function(e1,e2){
  return(e1@val<=e2)
})
setMethod("<=",c("numeric","fAD"),function(e1,e2){
  return(e1<=e2@val)
})
setMethod("<=",c("numeric","fAD2"),function(e1,e2){
  return(e1<=e2@val)
})






#
# binary numerical operations
#


setMethod("+",c("fAD","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val+e2@val,jac=e1@jac+e2@jac))
            }  
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })

setMethod("+",c("fAD2","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val+e2@val,jac=e1@jac+e2@jac,
                         hessian=e1@hessian+e2@hessian))
            }  
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })


setMethod("+",c("fAD","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val+e2,jac=e1@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })

setMethod("+",c("fAD2","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val+e2,jac=e1@jac,hessian=e1@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })

setMethod("+",c("numeric","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1+e2@val,jac=e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })

setMethod("+",c("numeric","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1+e2@val,jac=e2@jac,hessian=e2@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)+rep_len(e2,out.dim))
          })



setMethod("-",c("fAD","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val-e2@val,jac=e1@jac-e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })

setMethod("-",c("fAD2","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val-e2@val,jac=e1@jac-e2@jac,
                         hessian=e1@hessian-e2@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })

setMethod("-",c("fAD","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val-e2,jac=e1@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })

setMethod("-",c("fAD2","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val-e2,jac=e1@jac,hessian=e1@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })

setMethod("-",c("numeric","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1-e2@val,jac=-e2@jac))
            }
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })
setMethod("-",c("numeric","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1-e2@val,jac=-e2@jac,hessian=-e2@hessian))
            }
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)-rep_len(e2,out.dim))
          })



setMethod("*",c("fAD","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val*e2@val,jac=e2@val*e1@jac+e1@val*e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })

setMethod("*",c("fAD2","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              tmp <- e1@val*e2@hessian+e2@val*e1@hessian
              for(i in 1:length(e1)){
                tmp2 <- outer(e1@jac[i,],e2@jac[i,])
                tmp[i,,] <- tmp[i,,] + tmp2 + t(tmp2)
              }
              return(new("fAD2",val=e1@val*e2@val,jac=e2@val*e1@jac+e1@val*e2@jac,
                         hessian=tmp))
            }
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })



setMethod("*",c("fAD","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val*e2,jac=e2*e1@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })

setMethod("*",c("fAD2","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val*e2,jac=e2*e1@jac,hessian=e2*e1@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })


setMethod("*",c("numeric","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1*e2@val,jac=e1*e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })

setMethod("*",c("numeric","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1*e2@val,jac=e1*e2@jac,hessian=e1*e2@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)*rep_len(e2,out.dim))
          })


setMethod("/",c("fAD","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val/e2@val,jac=(1.0/e2@val)*e1@jac-(e1@val/(e2@val^2))*e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })

setMethod("/",c("fAD2","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              t1 <- 1.0/e2@val
              t2 <- e1@val/(e2@val^2)
              t3 <- 1.0/(e2@val^2)
              t4 <- 2.0*e1@val/(e2@val^3)
              tmp <- t1*e1@hessian - t2*e2@hessian
              for(i in 1:length(e1)){
                tmp2 <- outer(e1@jac[i,],e2@jac[i,])
                tmp[i,,] <- tmp[i,,] - t3[i]*(tmp2+t(tmp2)) + t4[i]*outer(e2@jac[i,],e2@jac[i,])  
              }
              return(new("fAD2",val=e1@val/e2@val,jac=t1*e1@jac-t2*e2@jac,hessian=tmp))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })



setMethod("/",c("fAD","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1@val/e2,jac=(1.0/e2)*e1@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })
setMethod("/",c("fAD2","numeric"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD2",val=e1@val/e2,jac=(1.0/e2)*e1@jac,hessian=(1.0/e2)*e1@hessian))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })


setMethod("/",c("numeric","fAD"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              return(new("fAD",val=e1/e2@val,jac=-(e1/(e2@val^2))*e2@jac))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })

setMethod("/",c("numeric","fAD2"),
          function(e1,e2){
            if(length(e1)==length(e2)){
              t1 <- (e1/(e2@val^2))
              t2 <- (2.0*e1/(e2@val^3))
              tmp <- -t1*e2@hessian
              for(i in 1:length(e1)) tmp[i,,] <- tmp[i,,] + t2[i]*outer(e2@jac[i,],e2@jac[i,])
              return(new("fAD2",val=e1/e2@val,jac=-t1*e2@jac,
                         hessian=tmp))
            } 
            out.dim <- max(length(e1),length(e2))
            return(rep_len(e1,out.dim)/rep_len(e2,out.dim))
          })




# lazy, consider better implementation
setMethod("^",c("fAD","fAD"),function(e1,e2){return(exp(e2*log(e1)))})
setMethod("^",c("fAD2","fAD2"),function(e1,e2){return(exp(e2*log(e1)))})
setMethod("^",c("fAD","numeric"),function(e1,e2){
  if(length(e1)==length(e2) || length(e2)==1){
    if(all(abs(e2-round(e2))<1.0e-14)){
      tmp <- e1@val^(e2-1)
      return(new("fAD",val=tmp*e1@val,jac=tmp*e2*e1@jac))
    } else {
      return(exp(e2*log(e1)))
    }
  } else {
    m <- max(length(e1),length(e2))
    return(rep_len(e1,m)^rep_len(e2,m))
  }
})
setMethod("^",c("fAD2","numeric"),function(e1,e2){
  
  if(length(e1)==length(e2) || length(e2)==1){
    
    if(all(abs(e2-round(e2))<1.0e-14)){ # e2 is integer
      
      tmp <- e1@val^(e2-1)
      t1 <- e1@val^(e2-2)*(e2^2-e2)
      tmp.h <- e2*tmp*e1@hessian
      for(i in 1:length(e1)){
        tmp.h[i,,] <- tmp.h[i,,] + t1[i]*outer(e1@jac[i,],e1@jac[i,])  
      }
      return(new("fAD2",val=tmp*e1@val,jac=tmp*e2*e1@jac,
                 hessian = tmp.h))
    } else {
      return(exp(e2*log(e1)))
    } 
  } else {
    m <- max(length(e1),length(e2))
    return(rep_len(e1,m)^rep_len(e2,m))
  }
})

setMethod("^",c("numeric","fAD"),function(e1,e2){return(exp(e2*log(e1)))})
setMethod("^",c("numeric","fAD2"),function(e1,e2){return(exp(e2*log(e1)))})

#
# Elementwise unary operations
#

setMethod("-",c("fAD"),function(e1){
  return(new("fAD",val=-e1@val,jac=-e1@jac))
})
setMethod("-",c("fAD2"),function(e1){
  return(new("fAD2",val=-e1@val,jac=-e1@jac,hessian=-e1@hessian))
})
setMethod("+",c("fAD"),function(e1){
  return(e1)
})
setMethod("+",c("fAD2"),function(e1){
  return(e1)
})


setMethod("exp",c("fAD"),
          function(x){
            tmp <- exp(x@val)
            return(new("fAD",val=tmp,jac=tmp*x@jac))
          })
setMethod("exp",c("fAD2"),
          function(x){
            tmp <- exp(x@val)
            tmp.h <- tmp*x@hessian
            for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] + tmp[i]*outer(x@jac[i,],x@jac[i,])
            return(new("fAD2",val=tmp,jac=tmp*x@jac,
                       hessian=tmp.h))
          })

setMethod("log",c("fAD"),
          function(x){
            return(new("fAD",val=log(x@val),jac=(1.0/x@val)*x@jac))
          })
setMethod("log",c("fAD2"),
          function(x){
            tmp <- log(x@val)
            t1 <- 1.0/x@val
            t2 <- -t1^2
            tmp.h <- t1*x@hessian
            for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] + t2[i]*outer(x@jac[i,],x@jac[i,])
            return(new("fAD2",val=tmp,jac=t1*x@jac,
                       hessian=tmp.h))
          })



setMethod("sqrt",c("fAD"),
          function(x){
            tmp <- sqrt(x@val)
            return(new("fAD",val=tmp,jac=(0.5/tmp)*x@jac))
          })
setMethod("sqrt",c("fAD2"),
          function(x){
            tmp <- sqrt(x@val)
            t1 <- 0.5/tmp
            t2 <- -0.5*t1/x@val
            tmp.h <- t1*x@hessian
            for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] + t2[i]*outer(x@jac[i,],x@jac[i,])
            return(new("fAD2",val=tmp,jac=t1*x@jac,
                       hessian=tmp.h))
          })



setMethod("lgamma",c("fAD"),
          function(x){
            return(new("fAD",val=lgamma(x@val),jac=digamma(x@val)*x@jac))
          })
setMethod("lgamma",c("fAD2"),
          function(x){
            tmp <- digamma(x@val)
            t1 <- trigamma(x@val)
            tmp.h <- tmp*x@hessian
            for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,]+t1[i]*outer(x@jac[i,],x@jac[i,])
            return(new("fAD2",val=lgamma(x@val),jac=tmp*x@jac,
                       hessian=tmp.h))
          })

#' @export
setGeneric("square",function(x) standardGeneric("square"))

#' @export
setMethod("square",c("numeric"),function(x){
  return(x^2)
})

#' @export
setMethod("square",c("fAD"),function(x){
  return(new("fAD",val=x@val^2,jac=2.0*x@val*x@jac))
})
#' @export
setMethod("square",c("fAD2"),function(x){
  t1 <- 2.0*x@val
  tmp.h <- t1*x@hessian
  for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] + 2.0*outer(x@jac[i,],x@jac[i,])
  return(new("fAD2",val=x@val^2,jac=t1*x@jac,
             hessian=tmp.h))
})


setMethod("sin",c("fAD"),function(x){
  new("fAD",val=sin(x@val),jac=cos(x@val)*x@jac)
})
setMethod("sin",c("fAD2"),function(x){
  t1 <- sin(x@val)
  t2 <- cos(x@val)
  tmp.h <- t2*x@hessian
  for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] - t1*outer(x@jac[i,],x@jac[i,])
  return(new("fAD2",val=t1,jac=t2*x@jac,hessian=tmp.h))
})

setMethod("cos",c("fAD"),function(x){
  new("fAD",val=cos(x@val),jac=-sin(x@val)*x@jac)
})
setMethod("cos",c("fAD2"),function(x){
  t.sin <- sin(x@val)
  t.cos <- cos(x@val)
  tmp.h <- -t.sin*x@hessian
  for(i in 1:length(x)) tmp.h[i,,] <- tmp.h[i,,] - t.cos*outer(x@jac[i,],x@jac[i,])
  return(new("fAD2",val=t.cos,jac=-t.sin*x@jac,hessian=tmp.h))
})


.AD_abs_core <- function(x){
  return((-1+2*as.numeric(x@val>=0))*x)
}
setMethod("abs",c("fAD"),.AD_abs_core)
setMethod("abs",c("fAD2"),.AD_abs_core)

#
# reductions etc
#

setMethod("sum",c("fAD"),function(x){
  return(new("fAD",val=sum(x@val),jac=matrix(colSums(x@jac),nrow=1)))
})
setMethod("sum",c("fAD2"),function(x){
  return(new("fAD2",val=sum(x@val),jac=matrix(colSums(x@jac),nrow=1),
             hessian=array(colSums(x@hessian),c(1,.GlobalEnv$fAD.nvars__,.GlobalEnv$fAD.nvars__))))
})
.AD_prod_core <- function(x,...,na.rm){
  tmp <- x[1]
  if(length(x)>1) for(i in 2:length(x)) tmp <- tmp*x[i]
  return(tmp)
}
setMethod("prod",c("fAD","ANY"),.AD_prod_core)
setMethod("prod",c("fAD2","ANY"),.AD_prod_core)







