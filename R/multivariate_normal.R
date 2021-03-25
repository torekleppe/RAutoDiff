

#' @export
mvnorm_lpdf_var <- function(x,mean,sigma=diag(length(mean))){
  if(!is.null(nrow(x))){
    xx <- t(x)
    nr <- nrow(xx)
    m <- TRUE
  } else {
    xx <- x
    nr <- length(xx)
    m <- FALSE
  }
  d <- length(mean)
  if(nr!=d) stop("dimension mismatch")
  dev <- xx-mean
  L <- cholL(sigma)
  SigInvDev <- backsolve(t(L),forwardsolve(L,dev))
  if(m) {
    biLinForm <- colSums(dev*SigInvDev)
  } else {
    biLinForm <- sum(dev*SigInvDev)
  }
  return(-d*0.918938533204673 - sum(log(diag(L))) - 0.5*biLinForm)
}

#' @export
mvnorm_lpdf_var_chol <- function(x,mean,L=diag(length(mean))){
  if(!is.null(nrow(x))){
    xx <- t(x)
    nr <- nrow(xx)
    m <- TRUE
  } else {
    xx <- x
    nr <- length(xx)
    m <- FALSE
  }
  d <- length(mean)
  if(nr!=d) stop("dimension mismatch")
  dev <- xx-mean
  SigInvDev <- backsolve(t(L),forwardsolve(L,dev))
  if(m) {
    biLinForm <- colSums(dev*SigInvDev)
  } else {
    biLinForm <- sum(dev*SigInvDev)
  }
  return(-d*0.918938533204673 - sum(log(diag(L))) - 0.5*biLinForm)
}

#' @export
mvnorm_lpdf_prec <- function(x,mean,Prec=diag(length(mean))){
  if(!is.null(nrow(x))){
    xx <- t(x)
    nr <- nrow(xx)
    m <- TRUE
  } else {
    xx <- x
    nr <- length(xx)
    m <- FALSE
  }
  d <- length(mean)
  if(nr!=d) stop("dimension mismatch")
  dev <- xx-mean

  SigInvDev <- Prec%*%dev
  
  if(m) {
    biLinForm <- colSums(dev*SigInvDev)
  } else {
    biLinForm <- sum(dev*SigInvDev)
  }
  return(-d*0.918938533204673 + sum(log(diag(cholL(Prec)))) - 0.5*biLinForm)
}

#' @export
mvnorm_lpdf_prec_chol <- function(x,mean,PL=diag(length(mean))){
  d <- length(mean)
  if(!is.null(nrow(x))){
    xx <- t(x)
    nr <- nrow(xx)
    m <- TRUE
  } else {
    xx <- x
    nr <- length(xx)
    m <- FALSE
  }
  
  if(nr!=d) stop("dimension mismatch")
  
  dev <- xx-mean
  
  tmp <- t(PL)%*%dev
  if(m) {
    biLinForm <- colSums(square(tmp))
  } else {
    biLinForm <- sum(square(tmp))
  }
  return(-d*0.918938533204673 + sum(log(diag(PL))) - 0.5*biLinForm)
}



