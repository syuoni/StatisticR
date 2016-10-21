#' Duration Models Estimation with MLE
#' 
#' @name mle_duration_models
#' @param t explained variable, vector, duration time
#' @param X explanatory variable in matrix
#' @param d explanatory variable, vector, whether the sample is failed (1=failure)
# 
#' @return a list with 5 elements, mle result
#'   \item{method}{'mle'}
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table}
NULL

exp.lnlike <- function(beta, args){
  # d is a vector (dim=n), each element of which indicates the sample is censored 
  # d=0 means censored
  t <- args$t
  d <- args$d
  X <- args$X
  
  # A %*% B always return a matrix, transfrom it to a vector
  Xb <- as.vector(X %*% beta)
  lnlike <- d*Xb - exp(Xb)*t
  return(-sum(lnlike))
}

exp.gr <- function(beta, args){
  t <- args$t
  d <- args$d
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  gr <- apply((d-exp(Xb)*t) * X, 2, sum)
  return(-gr)
}

#' @export 
#' @rdname mle_duration_models
mle.exp.estimate <- function(t, X, d){
  args <- list(t=t, X=as.matrix(X), d=d)
  args <- clean4regression(args)
  
  # with specifying gradient function, the result would be accurate (consistent to Stata)
  # without gradient function, a finite-difference approximation will be used
  # the initial params0 should not make likelihood function return a infinite value
  params0 <- rep(1e-5, dim(args$X)[2])
  names(params0) <- colnames(args$X)
  model.res <- mle.model(exp.lnlike, args, params0=params0, gr=exp.gr)
  return(model.res)
}

weibull.lnlike <- function(theta, args){
  # the last element in theta is lnp (scalar), 
  # and the others is beta (vector)
  n.theta <- length(theta)
  beta <- theta[1:n.theta-1]
  lnp <- theta[n.theta]
  
  # d is a vector (dim=n), each element of which indicates the sample is censored 
  # d=0 means censored
  t <- args$t
  d <- args$d
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  lnlike <- d*(Xb+lnp+(exp(lnp)-1)*log(t)) - exp(Xb)*(t**exp(lnp))
  return(-sum(lnlike))
}

weibull.gr <- function(theta, args){
  n.theta <- length(theta)
  beta <- theta[1:(n.theta-1)]
  lnp <- theta[n.theta]
  
  t <- args$t
  d <- args$d
  X <- args$X
  
  p <- exp(lnp)
  Xb <- as.vector(X %*% beta)
  beta.gr <- apply((d-exp(Xb)*(t**p)) * X, 2, sum)
  lnp.gr <- sum(d*(1+log(t)*p) - exp(Xb)*(t**p)*log(t)*p)
  return(-c(beta.gr, lnp.gr))
}

#' @export 
#' @rdname mle_duration_models
mle.weibull.estimate <- function(t, X, d){
  args <- list(t=t, X=as.matrix(X), d=d)
  args <- clean4regression(args)
  
  params0 <- rep(1e-5, dim(args$X)[2]+1)
  names(params0) <- c(colnames(args$X), 'lnp')
  model.res <- mle.model(weibull.lnlike, args, params0=params0, gr=weibull.gr)
  return(model.res)
}

gompertz.lnlike <- function(theta, args){
  # the last element in theta is gamma (scalar), 
  # and the others is beta (vector)
  n.theta <- length(theta)
  beta <- theta[1:n.theta-1]
  gamma <- theta[n.theta]
  
  t <- args$t
  d <- args$d
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  lnlike <- d*(Xb+gamma*t) - exp(Xb)*(exp(gamma*t)-1)/gamma
  return(-sum(lnlike))
}

gompertz.gr <- function(theta, args){
  n.theta <- length(theta)
  beta <- theta[1:n.theta-1]
  gamma <- theta[n.theta]
  
  t <- args$t
  d <- args$d
  X <- args$X
  
  Xb <- as.vector(X %*% beta)
  beta.gr <- apply((d-exp(Xb)*(exp(gamma*t)-1)/gamma) * X, 2, sum)
  gamma.gr <- sum(d*t - exp(Xb)*(exp(gamma*t)*(gamma*t-1)+1)/gamma**2)
  return(-c(beta.gr, gamma.gr))
}

#' @export 
#' @rdname mle_duration_models
mle.gompertz.estimate <- function(t, X, d){
  args <- list(t=t, X=as.matrix(X), d=d)
  args <- clean4regression(args)
  
  params0 <- rep(1e-5, dim(args$X)[2]+1)
  names(params0) <- c(colnames(args$X), 'gamma')
  model.res <- mle.model(gompertz.lnlike, args, params0=params0, gr=gompertz.gr)
  return(model.res)
}

