#' @title Linear Models Estimation with MLE
#' @author syuoni
#' 
#' @name mle_linear_models
#' @description Linear Models Estimation with MLE
#' @param y explained variable
#' @param X explanatory variable in matrix
#' 
#' @return a list with 4 elements, mle result
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table}
#' 
#' @export 

mle.linear.estimate <- function(y, X){
  args <- list(y=y, X=as.matrix(X))
  args <- clean4regression(args)
  
  # converge only when the initial params are specified accurately
  # BFGS method can converge, with warnings? but BFGS fails in Python
  # With specifying sigma <- abs(sigma), there would be no warnings
  # params0 <- c(1, 1, 1, -0.05, 0.5, 0.27)
  params0 <- rep(1, dim(args$X)[2]+1)
  names(params0) <- c(colnames(args$X), 'sigma')
  model.res <- mle.model(linear.lnlike, args, params0=params0)
  return(model.res)
}

linear.lnlike <- function(theta, args){
  n.theta <- length(theta)
  # 1:n-1 is (1:n)-1, namely c(0, 1, 2, ..., n-1)
  # 1:(n-1) is c(1, 2, 3, ..., n-1)
  beta <- theta[1:(n.theta-1)]
  sigma <- theta[n.theta]
  
  y <- args$y
  X <- args$X
  e <- y - X %*% beta
  
  # With specifying sigma <- abs(sigma), there would be no warnings
  # sigma may be estimated to be negative
  # general solution: sigma should be transformed to ln(sigma)
  sigma <- abs(sigma)
  ln.prob.density <- dnorm(e, mean=0., sd=sigma, log=TRUE)
  return(-sum(ln.prob.density))
}
