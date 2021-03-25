# note this function is limited to cases where the first
# argument is of type fAD
# any help on figuring out how to generalize most welcome

#' @export
setMethod("c",signature(x="ADtype"),function(x,...){
            # count length of out-vector
            d <- length(x)
            for(i in 1:(nargs()-1)){
              d <- d+length(...elt(i))
            }
            # allocate out-vector
            out <- .ADtype.zeros(d,class(x)[1]) # fAD.zeros(d)
            # fill out-vector
            n <- length(x)
            coeffRef(out,1:n) <- x
            
            for(i in 1:(nargs()-1)){
              o <- ...elt(i)
              n.o <- length(o)
              inds <- (n+1):(n+n.o)
              if(any(class(o)=="matrix")){
                out[inds] <- o[1:n.o]
              }
              if(any(class(o)=="numeric")){
                out[inds] <- o
              }
              if(any(class(o)==class(x)[1])){
                out[inds] <- o
              }
              if(any(class(o)=="AD_matrix")){
                out[inds] <- o@vals
              }
              n <- n+n.o
            }
            return(out)
          }
          )
