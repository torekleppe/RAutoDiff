
setGeneric("cholL",function(x) standardGeneric("cholL"))

#' @export
setMethod("cholL",c("matrix"),function(x){return(t(chol(x)))})

#' @export
setMethod("cholL",c("AD_matrix"),function(x){
  if(x@nrow != x@ncol) stop("must be square")
  if(! isSymmetric(x)) stop("must be symmetric")
  n <- x@nrow
  L <- .AD_matrix.zeros(n,n,class(x@vals)[1]) # fAD_matrix.zeros(n,n)
  for(i in 1:n){
    di <- i+n*(i-1)
    dd <- x@vals[di] - sum(square(.row.block(L,i)))
    if(dd@val[1]<=1.0e-14) stop("matrix not PD")
    L@vals[di] <- sqrt(dd)
    if(i<n){
      for(j in (i+1):n){
        odi <- j+n*(i-1)
        L@vals[odi] <- (x@vals[odi] 
                        - sum(.row.block(L,i)*.row.block(L,j)))/L@vals[di]
      }
    }
  }
  return(L) 
})

#' overload of the built in (upper-triangular) Cholesky 
#' @export
setMethod("chol",c("AD_matrix"),function(x,...){return(t(cholL(x)))})

setGeneric("solve.chol",function(a,b) standardGeneric("solve.chol"))


#' @export
setMethod("solve.chol",c("AD_matrix","ANY"),
          function(a,b){
            if(missing(b)) b <- diag(nrow(a))
            L <- cholL(a)
            t1 <- forwardsolve(L,b)
            return(backsolve(t(L),t1))
          })

#' @export
setMethod("solve.chol",c("matrix","ANY"),
          function(a,b){
            if(missing(b)) b <- diag(nrow(a))
            L <- cholL(a)
            return(backsolve(t(L),forwardsolve(L,b)))
          })

