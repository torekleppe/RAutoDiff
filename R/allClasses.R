
#' class representing a numeric vector and it's Jacobian matrix
setClass("fAD",
         slots = list(val="numeric",jac="matrix")
)
#' class representing a numeric vector, it's Jacobian matrix and an 3d array of Hessian matrices
setClass("fAD2",
         slots = list(val="numeric",jac="matrix",hessian="array")
)

.AD_types <- c("fAD","fAD2")

#' Union of available AD types
setClassUnion("ADtype",members=c("fAD","fAD2")
)

#' matrix class populated by AD types
setClass("AD_matrix",
         slots=list(nrow="integer",ncol="integer",vals="ADtype")
)


