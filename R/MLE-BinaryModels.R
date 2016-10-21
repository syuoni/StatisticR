#' Binary Models Estimation with MLE
#' 
#' @name mle_binary_models
#' @param y explained variable, binary
#' @param X explanatory variable in matrix
#' 
#' @return a list with 5 elements, mle result
#'   \item{method}{'mle'}
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table}
NULL

probit.lnlike <- function(beta, args){
  y <- args$y
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  lnlike <- y*pnorm(Xb, log.p=TRUE) + (1-y)*pnorm(-Xb, log.p=TRUE)
  return(-sum(lnlike))
}

probit.gr <- function(beta, args){
  y <- args$y
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  gr <- apply(y * dnorm(Xb)/pnorm(Xb)*X + (1-y) * dnorm(-Xb)/pnorm(-Xb)*(-X), 2, sum)
  return(-gr)
}

#' @export 
#' @rdname mle_binary_models
mle.probit.estimate <- function(y, X){
  args <- list(y=y, X=as.matrix(X))
  args <- clean4regression(args)
  
  params0 <- rep(1e-5, dim(args$X)[2])
  names(params0) <- colnames(args$X)
  model.res <- mle.model(probit.lnlike, args, params0=params0, gr=probit.gr)
  return(model.res)
}

logit.lnlike <- function(beta, args){
  y <- args$y
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  expXb = exp(Xb)
  lnlike <- y*log(expXb/(1+expXb)) + (1-y)*log(1/(1+expXb))
  return(-sum(lnlike))
}

logit.gr <- function(beta, args){
  y <- args$y
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  expXb = exp(Xb)
  gr <- apply(y * 1/(expXb+1)*X + (1-y) * (-expXb)/(expXb+1)*X, 2, sum)
  return(-gr)
}

#' @export 
#' @rdname mle_binary_models
mle.logit.estimate <- function(y, X){
  args <- list(y=y, X=as.matrix(X))
  args <- clean4regression(args)
  
  params0 <- rep(1e-5, dim(args$X)[2])
  names(params0) <- colnames(args$X)
  model.res <- mle.model(logit.lnlike, args, params0=params0, gr=logit.gr)
  return(model.res)
}
