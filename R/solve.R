
# implementation taken from Wikipedia, no pivoting done, so it should only
# be used on small, well-scaled problems
.LU_crout_core <- function(A){
  n <- nrow(A)
  if(n!=ncol(A)) stop("matrix must be square")
  L <- .AD_matrix.zeros(n,n,class(A@vals)[1])
  U <- .AD_matrix.zeros(n,n,class(A@vals)[1])
  for(i in 1:n){
    L[i,1] <- A[i,1]
    U[i,i] <- 1.0
  }
  U[1,2:n] <- .row.block(A,1,2,n)/L[1,1]
  for(i in 2:n){
    for(j in 2:i){
      L[i,j] <- A[i,j] - .row.block(L,i,1,j-1)%*%.col.block(U,j,1,j-1)
    }
    if(i<n){
      for(j in (i+1):n){
        U[i,j] <- (A[i,j] - .row.block(L,i,1,i-1) %*% .col.block(U,j,1,i-1))/L[i,i]
      }
    }
  }
  return(list(L=L,U=U))
}


# #' @export
#LUdecomp <- function(A){return(.LU_crout_core(A))}


#' Overload of the solve()-function for AD-types
#' @export
setMethod("solve",c("AD_matrix","ANY"),
          function(a,b,...){
            if(missing(b)) b <- diag(nrow(a))
            lu <- .LU_crout_core(a)
            t1 <- forwardsolve(lu$L,b)
            tmp <- backsolve(lu$U,t1)
            
            
            if(class(tmp)[1]=="AD_matrix"){
              if(ncol(tmp)==1){
                return(tmp@vals)
              } 
            }
            return(tmp)
              
          })
#' Overload of the solve()-function for AD-types
#' @export
setMethod("solve",c("matrix","ADtype"),
          function(a,b,...){
            return((solve(a,...)%*%b)@vals)
            })
#' Overload of the solve()-function for AD-types
#' @export
setMethod("solve",c("matrix","AD_matrix"),
          function(a,b,...){
            return(solve(a,...)%*%b)
          })

#' Deterimant of AD_matrix type
#' @export
setMethod("determinant",c("AD_matrix","ANY"),
          function(x,logarithm=TRUE,...){
            lu <- .LU_crout_core(x)
            if(logarithm){
              return(sum(log(diag(lu$L))))
            } else {
              return(prod(diag(lu$L)))
            }
          })

#' Deterimant of AD_matrix type
#' @export
setMethod("det",c("AD_matrix"),function(x,...){return(determinant(x,logarithm = FALSE))})