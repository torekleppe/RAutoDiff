

setMethod("show",c("AD_matrix"),
          function(object){
            print(paste0(object@nrow, " by ",object@ncol, " fAD_matrix, values : "),quote = FALSE)
            print(matrix(object@vals@val,nrow=object@nrow),quote=FALSE)
            print(paste0("with ",class(object@vals)[1]," AD types"),quote = FALSE)
          })

#' @export
setMethod("nrow",c("AD_matrix"),function(x){return(x@nrow)})

#' @export
setMethod("ncol",c("AD_matrix"),function(x){return(x@ncol)})

#' @export
setMethod("dim",c("AD_matrix"),function(x){return(c(x@nrow,x@ncol))})

.AD_matrix.zeros <- function(nrow=1,ncol=1,type){
  return(new("AD_matrix",nrow=as.integer(nrow),ncol=as.integer(ncol),
             vals=.ADtype.zeros(nrow*ncol,type)))
}

# #' @export
#fAD_matrix.zeros <- function(nrow=1,ncol=1){
#  return(new("AD_matrix",nrow=as.integer(nrow),ncol=as.integer(ncol),vals=fAD.zeros(nrow*ncol)))
#}

# #' @export
#fAD_matrix.identity <- function(n=1){
#  return(new("AD_matrix",nrow=as.integer(n),ncol=as.integer(n),vals=.as.fAD(diag(n)[1:(n^2)])))
#}


#' @export
setMethod("matrix",c("ADtype","ANY","ANY"),
          function (data, nrow=NULL , ncol=NULL) {
            
            
            if(is.null(nrow) && is.null(ncol)){
              nrow <- length(data)
              ncol <- 1
            } else if(is.null(nrow) && length(data)>1){
              nrow <- length(data)/ncol
              if(nrow-floor(nrow)>1.0e-10) stop("dimension mismatch")
            } else if(is.null(ncol) && length(data)>1){
              ncol <- length(data)/nrow
              if(ncol-floor(ncol)>1.0e-10) stop("dimension mismatch")
            }
            
            nrow <- as.integer(nrow)
            ncol <- as.integer(ncol)
            
            if(length(data)==nrow*ncol){
              return(new("AD_matrix",nrow=nrow,ncol=ncol,vals=data))
            } else if(length(data)==1){
              return(new("AD_matrix",nrow=nrow,ncol=ncol,vals=rep.int(data,nrow*ncol)))
            } else if(missing(nrow) && missing(ncol)){
              return(new("AD_matrix",nrow=length(data),ncol=1L,vals=data))
            } else {
              stop("dimension mismatch")
            }
          })

setMethod("[",c("AD_matrix","ANY","ANY","missing"),
          function (x, i, j, ..., drop = TRUE){
            
            if(nargs()==2) {
              if(! missing(i)){
                return(x@vals[i])
              } else {
                return(x)
              }
            } else if(nargs()==3){
              if(missing(i)) i <- 1:x@nrow
              if(missing(j)) j <- 1:x@ncol
              nr <- x@nrow
              inds <- as.vector(outer(i,j,function(x,y){return(x+(y-1)*nr)}))
              if(length(i)==1 || length(j)==1) return(x@vals[inds])
              return(new("AD_matrix",nrow=length(i),ncol=length(j),vals=x@vals[inds]))
            }
            stop()
          })

setMethod("[<-",signature(x = "AD_matrix","ANY","ANY","ANY"),
          function (x, i, j, ..., value){
            if(nargs()==3){
              if(missing(i)) i <- 1:(x@nrow*x@ncol)
              if(length(i) != length(value)) stop("dimension mismatch") 
              x@vals[i] <- value
              return(x)
            } else {
              if(missing(i)) i <- 1:(x@nrow)
              if(missing(j)) j <- 1:(x@ncol)
              nr <- x@nrow
              inds <- as.vector(outer(i,j,function(x,y){return(x+(y-1)*nr)}))
              x@vals[inds] <- value
              return(x)  
            }
          })

#' @export
setGeneric(".col.block",function(x,which.col,from=1,to=x@nrow) 
  standardGeneric(".col.block"))

setMethod(".col.block",c(x="AD_matrix","ANY","ANY","ANY"),
          function(x,which.col,from=1,to=x@nrow){
            col.start <- x@nrow*(which.col-1)
            return(x@vals[(col.start+from):(col.start+to)])
          })

#' @export
setGeneric(".row.block",function(x,which.row,from=1,to=x@nrow) 
  standardGeneric(".row.block"))

setMethod(".row.block",c(x="AD_matrix","ANY","ANY","ANY"),
          function(x,which.row,from=1,to=x@ncol){
            return(x@vals[which.row+((from-1):(to-1))*x@nrow])
          })
# mainly for test purposes
#' @export
setGeneric(".as.double.matrix",function(x) standardGeneric(".as.double.matrix"))
setMethod(".as.double.matrix",c("AD_matrix"),function(x){
  return(matrix(x@vals@val,x@nrow,x@ncol))
})


#' @export
setMethod("length",c("AD_matrix"),function(x){return(length(x@vals))})


#
# Element-wise binary operations
#

setMethod("+",c("AD_matrix","AD_matrix"),function(e1,e2) {
  if(e1@nrow!=e2@nrow || e1@ncol!=e2@ncol) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals+e2@vals))
})
setMethod("+",c("AD_matrix","matrix"),function(e1,e2) {
  if(e1@nrow!=dim(e2)[1] || e1@ncol!=dim(e2)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals+e2[1:length(e2)]))
})
setMethod("+",c("matrix","AD_matrix"),function(e1,e2) {
  if(e2@nrow!=dim(e1)[1] || e2@ncol!=dim(e1)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,vals=e1[1:length(e1)]+e2@vals))
})
setMethod("+",c("AD_matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals+rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("+",c("ADtype","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))+e2@vals))
})
setMethod("+",c("AD_matrix","numeric"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals+rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("+",c("numeric","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))+e2@vals))
})


setMethod("-",c("AD_matrix","AD_matrix"),function(e1,e2) {
  if(e1@nrow!=e2@nrow || e1@ncol!=e2@ncol) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals-e2@vals))
})
setMethod("-",c("AD_matrix","matrix"),function(e1,e2) {
  if(e1@nrow!=dim(e2)[1] || e1@ncol!=dim(e2)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals-e2[1:length(e2)]))
})
setMethod("-",c("matrix","AD_matrix"),function(e1,e2) {
  if(e2@nrow!=dim(e1)[1] || e2@ncol!=dim(e1)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,vals=e1[1:length(e1)]-e2@vals))
})
setMethod("-",c("AD_matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals-rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("-",c("ADtype","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))-e2@vals))
})
setMethod("-",c("AD_matrix","numeric"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals-rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("-",c("numeric","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))-e2@vals))
})

setMethod("*",c("AD_matrix","AD_matrix"),function(e1,e2) {
  if(e1@nrow!=e2@nrow || e1@ncol!=e2@ncol) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals*e2@vals))
})
setMethod("*",c("AD_matrix","matrix"),function(e1,e2) {
  if(e1@nrow!=dim(e2)[1] || e1@ncol!=dim(e2)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals*e2[1:length(e2)]))
})
setMethod("*",c("matrix","AD_matrix"),function(e1,e2) {
  if(e2@nrow!=dim(e1)[1] || e2@ncol!=dim(e1)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,vals=e1[1:length(e1)]*e2@vals))
})
setMethod("*",c("AD_matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals*rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("*",c("ADtype","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))*e2@vals))
})
setMethod("*",c("AD_matrix","numeric"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals*rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("*",c("numeric","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))*e2@vals))
})

setMethod("/",c("AD_matrix","AD_matrix"),function(e1,e2) {
  if(e1@nrow!=e2@nrow || e1@ncol!=e2@ncol) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals/e2@vals))
})
setMethod("/",c("AD_matrix","matrix"),function(e1,e2) {
  if(e1@nrow!=dim(e2)[1] || e1@ncol!=dim(e2)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals/e2[1:length(e2)]))
})
setMethod("/",c("matrix","AD_matrix"),function(e1,e2) {
  if(e2@nrow!=dim(e1)[1] || e2@ncol!=dim(e1)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,vals=e1[1:length(e1)]/e2@vals))
})
setMethod("/",c("AD_matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals/rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("/",c("ADtype","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))/e2@vals))
})
setMethod("/",c("AD_matrix","numeric"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals/rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("/",c("numeric","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))/e2@vals))
})


setMethod("^",c("AD_matrix","AD_matrix"),function(e1,e2) {
  if(e1@nrow!=e2@nrow || e1@ncol!=e2@ncol) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals^e2@vals))
})
setMethod("^",c("AD_matrix","matrix"),function(e1,e2) {
  if(e1@nrow!=dim(e2)[1] || e1@ncol!=dim(e2)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=e1@vals^e2[1:length(e2)]))
})
setMethod("^",c("matrix","AD_matrix"),function(e1,e2) {
  if(e2@nrow!=dim(e1)[1] || e2@ncol!=dim(e1)[2]) stop("dimension mismatch")
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,vals=e1[1:length(e1)]^e2@vals))
})
setMethod("^",c("AD_matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals^rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("^",c("ADtype","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))^e2@vals))
})
setMethod("^",c("AD_matrix","numeric"),function(e1,e2){
  return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,
             vals=e1@vals^rep_len(e2,length.out=(e1@nrow*e1@ncol))))
})
setMethod("^",c("numeric","AD_matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=e2@nrow,ncol=e2@ncol,
             vals=rep_len(e1,length.out=(e2@nrow*e2@ncol))^e2@vals))
})


#
# element-wise expressions of the type fAD op matrix and matrix op fAD, both resulting in
# fAD_matrix types
#
setMethod("+",c("ADtype","matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e2),ncol=ncol(e2),vals=e1+e2[1:length(e2)]))
})
setMethod("+",c("matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e1),ncol=ncol(e1),vals=e1[1:length(e1)]+e2))
})
setMethod("-",c("ADtype","matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e2),ncol=ncol(e2),vals=e1-e2[1:length(e2)]))
})
setMethod("-",c("matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e1),ncol=ncol(e1),vals=e1[1:length(e1)]-e2))
})
setMethod("*",c("ADtype","matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e2),ncol=ncol(e2),vals=e1*e2[1:length(e2)]))
})
setMethod("*",c("matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e1),ncol=ncol(e1),vals=e1[1:length(e1)]*e2))
})
setMethod("/",c("ADtype","matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e2),ncol=ncol(e2),vals=e1/e2[1:length(e2)]))
})
setMethod("/",c("matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e1),ncol=ncol(e1),vals=e1[1:length(e1)]/e2))
})
setMethod("^",c("ADtype","matrix"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e2),ncol=ncol(e2),vals=e1^e2[1:length(e2)]))
})
setMethod("^",c("matrix","ADtype"),function(e1,e2){
  return(new("AD_matrix",nrow=nrow(e1),ncol=ncol(e1),vals=e1[1:length(e1)]^e2))
})

#
# Element-wise unary operations
#
setMethod("-",c("AD_matrix"),
          function(e1){
            return(new("AD_matrix",nrow=e1@nrow,ncol=e1@ncol,vals=-e1@vals))
          })
setMethod("+",c("AD_matrix"),
          function(e1){
            return(e1)
          })

setMethod("exp",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=exp(x@vals)))
          })
setMethod("log",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=log(x@vals)))
          })
setMethod("sqrt",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=sqrt(x@vals)))
          })
setMethod("lgamma",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=lgamma(x@vals)))
          })
setMethod("sin",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=sin(x@vals)))
          })
setMethod("cos",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=cos(x@vals)))
          })
setMethod("abs",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=abs(x@vals)))
          })

#' @export
setMethod("square",c("AD_matrix"),
          function(x){
            return(new("AD_matrix",nrow=x@nrow,ncol=x@ncol,vals=square(x@vals)))
          })

#
# Basic linear algebra
#

# matrix-matrix products
setMethod("%*%",c("AD_matrix","AD_matrix"),
          function(x,y){
            if(x@ncol != y@nrow) stop("non-conformable arguments")
            out.nrow <- x@nrow
            out.ncol <- y@ncol
            ret <- .AD_matrix.zeros(out.nrow,out.ncol,class(x@vals)[1]) #fAD_matrix.zeros(out.nrow,out.ncol) 
            k <- 1
            for(j in 1:out.ncol){
              for(i in 1:out.nrow){
                ret@vals[k] <- sum(.row.block(x,i)*.col.block(y,j))
                k <- k+1
              }
            }
            return(ret)
          })

setMethod("%*%",c("AD_matrix","matrix"),
          function(x,y){
            if(x@ncol != dim(y)[1]) stop("non-conformable arguments")
            out.nrow <- x@nrow
            out.ncol <- dim(y)[2]
            ret <- .AD_matrix.zeros(out.nrow,out.ncol,class(x@vals)[1]) #fAD_matrix.zeros(out.nrow,out.ncol) 
            k <- 1
            for(j in 1:out.ncol){
              for(i in 1:out.nrow){
                ret@vals[k] <- sum(.row.block(x,i)*y[,j])
                k <- k+1
              }
            }
            return(ret)
          })

setMethod("%*%",c("matrix","AD_matrix"),
          function(x,y){
            if(dim(x)[2] != y@nrow) stop("non-conformable arguments")
            out.nrow <- dim(x)[1]
            out.ncol <- y@ncol
            ret <- .AD_matrix.zeros(out.nrow,out.ncol,class(y@vals)[1]) # fAD_matrix.zeros(out.nrow,out.ncol) 
            k <- 1
            for(j in 1:out.ncol){
              for(i in 1:out.nrow){
                ret@vals[k] <- sum(x[i,]*.col.block(y,j))
                k <- k+1
              }
            }
            return(ret)
          })
#
#matrix-vector products
#
setMethod("%*%",c("AD_matrix","ADtype"), # note; somewhat lazy
          function(x,y){
            return(x%*%new("AD_matrix",nrow=length(y),ncol=1L,vals=y))
          })
setMethod("%*%",c("ADtype","AD_matrix"),
          function(x,y){
            return(new("AD_matrix",nrow=1L,ncol=length(x),vals=x)%*%y)
          })
setMethod("%*%",c("AD_matrix","numeric"),
          function(x,y){
            return(x%*%matrix(y,length(y),1))
          })
setMethod("%*%",c("numeric","AD_matrix"),
          function(x,y){
            return(matrix(x,1,length(x))%*%y)
          })

setMethod("%*%",c("matrix","ADtype"),
          function(x,y){
            return(x%*%matrix(y,nrow=length(y),ncol=1))
          })
setMethod("%*%",c("ADtype","matrix"),
          function(x,y){
            return(matrix(x,nrow=1,ncol=length(x))%*%y)
          })
#
# dot-products
#
setMethod("%*%",c("ADtype","ADtype"),function(x,y){return(sum(x*y))})

setMethod("%*%",c("ADtype","numeric"),function(x,y){return(sum(x*y))})

setMethod("%*%",c("numeric","ADtype"),function(x,y){return(sum(x*y))})



# note, only check the value, not the AD-info

#' @export
setMethod("isSymmetric",c("AD_matrix"),function(object,...){
  return(isSymmetric(matrix(object@vals@val,nrow=object@nrow,ncol=object@ncol)))
})



# note explicit copy, possible room for improvement

#' @export
setMethod("t",c("AD_matrix"),function(x){
  
  ret <- .AD_matrix.zeros(x@ncol,x@nrow,class(x@vals)[1]) # fAD_matrix.zeros(x@ncol,x@nrow)
  rind <- 1:x@nrow
  for(j in 1:x@ncol){
    ret@vals[j+(x@ncol*(rind-1))] <- x@vals[rind+(x@nrow*(j-1))]
  }
  return(ret)
})

#' @export
setMethod("diag",
          signature(x = "AD_matrix","ANY","ANY","ANY"),
          function (x = 1, nrow, ncol, names = TRUE) 
          {
            n <- min(x@nrow,x@ncol)
            return(x@vals[(1:n) + (0:(n-1))*x@nrow])
          })

#' @export
setMethod("diag",
          signature(x = "ADtype","ANY","ANY","ANY"),
          function (x = 1, nrow, ncol, names = TRUE) 
          {
            n <- length(x)
            M <- .AD_matrix.zeros(n,n,class(x)[1]) # fAD_matrix.zeros(n,n)
            inds <- 1:n
            M@vals[inds+(inds-1)*n] <- x
            return(M)
          })


#
# forward/back-solve core methods
#
.AD_matrix.backsolve_core <- function(r,b){
  n <- length(b)
  x <- .ADtype.zeros(n,class(r@vals)[1],class(b)[1]) #fAD.zeros(n)
  coeffRef(x,n) <- coeff(b,n)/coeff(r@vals,n*n)
  for(i in (n-1):1){
    coeffRef(x,i) <- 
      (coeff(b,i)-sum(coeff(x,(i+1):n)*.row.block(r,i,i+1,n)))/coeff(r@vals,i+(i-1)*n)
  }
  return(x)
}

.matrix.backsolve_core <- function(r,b){
  n <- length(b)
  x <- .ADtype.zeros(n,class(b)[1]) #fAD.zeros(n)
  coeffRef(x,n) <- coeff(b,n)/r[n,n]
  for(i in (n-1):1){
    coeffRef(x,i) <- 
      (coeff(b,i)-sum(coeff(x,(i+1):n)*r[i,(i+1):n]))/r[i,i]
  }
  return(x)
}

.AD_matrix.forwardsolve_core <- function(l,b){
  n <- length(b)
  x <- .ADtype.zeros(n,class(l@vals)[1],class(b)[1]) #fAD.zeros(n)
  coeffRef(x,1) <- coeff(b,1)/coeff(l@vals,1)
  for(i in 2:n){
    coeffRef(x,i) <- (coeff(b,i) - sum(.row.block(l,i,1,i-1)*coeff(x,1:(i-1))))/coeff(l@vals,i+n*(i-1))
  }
  return(x)
}

.matrix.forwardsolve_core <- function(l,b){
  n <- length(b)
  x <- .ADtype.zeros(n,class(b)[1]) #  fAD.zeros(n)
  coeffRef(x,1) <- coeff(b,1)/l[1,1]
  for(i in 2:n){
    coeffRef(x,i) <- (coeff(b,i) - sum(l[i,1:(i-1)]*coeff(x,1:(i-1))))/l[i,i]
  }
  return(x)
}

#
#solves u%*%y=x for y when u is upper-triangular
#

#' @export
setMethod("backsolve",c("AD_matrix","AD_matrix"),
          function(r, x){
            nsol <- ncol(x)
            n <- nrow(r)
            y <- .AD_matrix.zeros(n,nsol,class(r@vals)[1]) # fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .AD_matrix.backsolve_core(r,.col.block(x,i))
            return(y)
          })

#' @export
setMethod("backsolve",c("AD_matrix","matrix"),
          function(r, x){
            nsol <- ncol(x)
            n <- nrow(r)
            y <- .AD_matrix.zeros(n,nsol,class(r@vals)[1]) # fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .AD_matrix.backsolve_core(r,x[,i])
            return(y)
          })

#' @export
setMethod("backsolve",c("AD_matrix","ADtype"),
          function(r,x){
            return(.AD_matrix.backsolve_core(r,x))
          })

#' @export
setMethod("backsolve",c("AD_matrix","numeric"),
          function(r,x){
            return(.AD_matrix.backsolve_core(r,x))
          })

#' @export
setMethod("backsolve",c("matrix","AD_matrix"),
          function(r, x){
            nsol <- ncol(x)
            n <- nrow(r)
            y <- .AD_matrix.zeros(n,nsol,class(x@vals)[1]) #fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .matrix.backsolve_core(r,.col.block(x,i))
            return(y)
          })

#' @export
setMethod("backsolve",c("matrix","ADtype"),
          function(r,x){
            return(.matrix.backsolve_core(r,x))
          })

#
# solves l%*%y = x for y when l is lower triangular
#

#' @export
setMethod("forwardsolve",c("AD_matrix","AD_matrix"),
          function(l,x){
            nsol <- ncol(x)
            n <- nrow(l)
            y <- .AD_matrix.zeros(n,nsol,class(l@vals)[1])  #fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .AD_matrix.forwardsolve_core(l,.col.block(x,i))
            return(y)
          })

#' @export
setMethod("forwardsolve",c("AD_matrix","matrix"),
          function(l,x){
            nsol <- ncol(x)
            n <- nrow(l)
            y <- .AD_matrix.zeros(n,nsol,class(l@vals)[1]) #fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .AD_matrix.forwardsolve_core(l,x[,i])
            return(y)
          })

#' @export
setMethod("forwardsolve",c("AD_matrix","ADtype"),
          function(l,x){
            return(.AD_matrix.forwardsolve_core(l,x))
          })

#' @export
setMethod("forwardsolve",c("AD_matrix","numeric"),
          function(l,x){
            return(.AD_matrix.forwardsolve_core(l,x))
          })

#' @export
setMethod("forwardsolve",c("matrix","AD_matrix"),
          function(l,x){
            nsol <- ncol(x)
            n <- nrow(l)
            y <- .AD_matrix.zeros(n,nsol,class(x@vals)[1]) #fAD_matrix.zeros(n,nsol)
            for(i in 1:nsol) coeffRef(y@vals,1:n+(i-1)*n) <- .matrix.forwardsolve_core(l,.col.block(x,i))
            return(y)
          })

#' @export
setMethod("forwardsolve",c("matrix","ADtype"),
          function(l,x){
            return(.matrix.forwardsolve_core(l,x))
          })




# general solve methods are based on explicitly calculating the
# inverse of the matrix, as this leads to simple derivative formulas
# .fAD_matrix.inverse_core <- function(x,...){
#   xi.d <- solve(.as.double.matrix(x),...)
#   n <- nrow(xi.d)
#   # values
#   vals <- fAD.zeros(n^2)
#   vals@val <- xi.d[1:length(xi.d)]
#   
#   # derivatives 
#   for(l in 1:n){
#     for(k in 1:n){
#       out.ind <- k + n*(l-1)
#       for(j in 1:n){
#         for(i in 1:n){
#           vals@jac[out.ind,] <- vals@jac[out.ind,] - xi.d[k,i]*xi.d[j,l]*x@vals@jac[i+n*(j-1),]
#         }
#       }
#     }
#   }
#   return(new("AD_matrix",nrow=n,ncol=n,vals=vals))
# }

#
# reductions
#

#' @export
setMethod("sum","AD_matrix",function(x){return(sum(x@vals))})

#' @export
setMethod("colSums","AD_matrix",function(x){
  n <- ncol(x)
  ret <- .ADtype.zeros(n,class(x@vals)[1])
  for(j in 1:n){
    coeffRef(ret,j) <- sum(.col.block(x,j))
  }
  return(ret)
})

#' @export
setMethod("rowSums","AD_matrix",function(x){
  n <- nrow(x)
  ret <- .ADtype.zeros(n,class(x@vals)[1])
  for(j in 1:n){
    coeffRef(ret,j) <- sum(.row.block(x,j))
  }
  return(ret)
})



