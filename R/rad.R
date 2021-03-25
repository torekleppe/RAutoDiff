
#source("RADv.R")
#source("RADm.R")
#source("combine.R")
#source("solve.R")
#source("chol.R")


#Todo: implement these with generics, in particular to enable
# handling of matrices



#' Set independent variable for first order derivatives
#'
#' @param x numeric vector to be differentiated with respect to
#' @return a fAD vector with value being x and jacobian being the identity
#' @examples 
#' x <- independent(c(0,1))
#' y <- sum(square(x))/2
#' print(gradient(y)) # should be c(0,1)
#' 
#' # define some function
#' f <- function(x){ 
#' M <- matrix(x,2,2)
#' return(sum(solve(M,c(-1,2))))
#' }
#' 
#' # run function with regular variables
#' print(f(c(1,-1,2,4)))
#' 
#' # run function with AD-variables
#' x <- independent(c(1,-1,2,4))
#' y <- f(x)
#' print(value(y))
#' print(gradient(y))
#' 
#' @export
independent <-function(x){
  .GlobalEnv$fAD.nvars__ <- length(x)
  return(new("fAD",val=x,jac=diag(length(x))))
}

#' Set independent variable for first and second order derivatives
#'
#' @param x numeric vector to be differentiated with respect to
#' @return a fAD2 vector with value being x, jacobian being the identity, and hessian being zero
#' @examples 
#' x <- independent2(c(0,1))
#' y <- sum(square(x))/2
#' print(gradient(y)) # should be c(0,1)
#' print(hessian(y)) # should be identity
#'  
#' # define some function
#' f <- function(x){ 
#' M <- matrix(x,2,2)
#' return(sum(solve(M,c(-1,2))))
#' }
#' 
#' # run function with regular variables
#' print(f(c(1,-1,2,4)))
#' 
#' # run function with second order AD-variables
#' x <- independent2(c(1,-1,2,4))
#' y <- f(x)
#' print(value(y))
#' print(gradient(y))
#' print(hessian(y))
#' 
#' @export
independent2 <- function(x){
  n <- length(x)
  .GlobalEnv$fAD.nvars__ <- n
  return(new("fAD2",val=x,jac=diag(n),hessian=array(0.0,c(n,n,n))))
}

#' Get value from AD type
#'
#' @param x An ADtype or AD_matrix type (overloaded also for numeric and matrix types)
#' @return The value of the AD type
#' @examples 
#' # define some function
#' f <- function(x){ 
#' M <- matrix(x,2,2)
#' return(sum(solve(M,c(-1,2))))
#' }
#' 
#' # run function with regular variables
#' print(f(c(1,-1,2,4)))
#' 
#' # run function with second order AD-variables
#' x <- independent2(c(1,-1,2,4))
#' y <- f(x)
#' print(value(y))
#' print(gradient(y))
#' print(hessian(y))
#' @export
setGeneric("value",function(y) standardGeneric("value") )

#' @export
setMethod("value","fAD",function(y){
  return(y@val)
})
#' @export
setMethod("value","fAD2",function(y){
  return(y@val)
})

#' @export
setMethod("value","AD_matrix",function(y){
  return(.as.double.matrix(y))
})

#' @export
setMethod("value","numeric",function(y){return(y)})

#' @export
setMethod("value","matrix",function(y){return(y)})

#' Get gradient from scalar AD type
#'
#' @param x An ADtype or AD_matrix type (overloaded also for numeric and matrix types)
#' @return The gradient of the (scalar) AD type
#' @examples 
#' # define some function
#' f <- function(x){ 
#' M <- matrix(x,2,2)
#' return(sum(solve(M,c(-1,2))))
#' }
#' 
#' # run function with regular variables
#' print(f(c(1,-1,2,4)))
#' 
#' # run function with second order AD-variables
#' x <- independent2(c(1,-1,2,4))
#' y <- f(x)
#' print(value(y))
#' print(gradient(y))
#' print(hessian(y))
#' @export
setGeneric("gradient",function(y) standardGeneric("gradient"))

#' @export
setMethod("gradient","fAD", function(y){
  if(length(y)>1) stop("gradient requires scalar argument")
  return(as.vector(y@jac))
})
#' @export
setMethod("gradient","fAD2", function(y){
  if(length(y)>1) stop("gradient requires scalar argument")
  return(as.vector(y@jac))
})
#' @export
setMethod("gradient","AD_matrix", function(y){
  if(nrow(y)>1 || ncol(y)>1) stop("gradient requires 1x1 matrix")
  return(as.vector(y@vals@jac))
})



setGeneric("jacobian",function(y) standardGeneric("jacobian"))

#' @export
setMethod("jacobian","fAD",function(y){
  return(y@jac)
})
#' @export
setMethod("jacobian","fAD2",function(y){
  return(y@jac)
})


setGeneric("hessian",function(y) standardGeneric("hessian"))

#' @export
setMethod("hessian","fAD2",function(y){
  if(length(y)>1) stop("hessian requires scalar argument")
  return(matrix(y@hessian,.GlobalEnv$fAD.nvars__))
})
#' @export
setMethod("hessian","AD_matrix",function(y){
  if(nrow(y)>1 || ncol(y)>1) stop("hessian requires 1x1 argument")
  if(class(y@vals)[1] != "fAD2") stop("hessian can only be extracted after call to independent2()")
  return(matrix(y@vals@hessian[1,,],.GlobalEnv$fAD.nvars__))
})



# fun <- function(f){
#   m <- matrix(c(f[1],f[2],f[2],f[3]),2,2)
#   f <- (c(f[2],2)%*%solve(m,c(-1,2)))
#   return(f)
#   
# }
# 
# x <- c(1.0,1.5,4.0)
# 
# 
# ax <- independent(x)
# g0 <- fun(x)
# fd.grad <- 0*x
# h <- 1.0e-5
# for(i in 1:length(x)){
#   xx <- x
#   xx[i] <- x[i]+h 
#   fd.grad[i] <- (fun(xx)-g0)/h
# }
# ag <- fun(ax)
# print(gradient(ag))
# print(fd.grad)
# 
# 
