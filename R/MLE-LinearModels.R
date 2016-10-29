#' Linear Models Estimation with MLE
#' 
#' @name mle_linear_models
#' @param y explained variable
#' @param X explanatory variable in matrix
#' 
#' @return a list with 5 elements, mle result
#'   \item{method}{'mle'}
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table} 
#' @export 

mle.linear.estimate <- function(y, X){
  args <- list(y=y, X=as.matrix(X))
  args <- clean4regression(args)
  
  # converge only when the initial params are specified accurately
  # BFGS method can converge, with warnings? but BFGS fails in Python
  # With specifying sigma <- abs(sigma), there would be no warnings
  # params0 <- c(1, 1, 1, -0.05, 0.5, 0.27)
  params0 <- rep(1, dim(args$X)[2]+1)
  names(params0) <- c(colnames(args$X), 'lnsigma')
  model.res <- mle.model(linear.lnlike, args, params0=params0, gr=linear.gr)
  return(model.res)
}

linear.lnlike <- function(theta, args){
  n.theta <- length(theta)
  # 1:n-1 is (1:n)-1, namely c(0, 1, 2, ..., n-1)
  # 1:(n-1) is c(1, 2, 3, ..., n-1)
  # Note that sigma > 0 and ln(sigma) has no restrict
  # if the MLE result for any parameter theta, is theta_hat, 
  # then the MLE result for f(theta) is f(theta_hat), namyly f(theta)_hat = f(theta_hat)
  beta <- theta[1:(n.theta-1)]
  lnsigma <- theta[n.theta]
  
  y <- args$y
  X <- args$X
  Xb <- as.vector(X %*% beta)
  
  # With specifying sigma > 0, there would be no warnings, else sigma may be estimated to be negative
  # general solution: sigma should be transformed to ln(sigma)
  sigma <- exp(lnsigma)
  lnlike <- -0.5*(y-Xb)**2/sigma**2 - 0.5*log(2*pi) - lnsigma
  return(-sum(lnlike))
}

linear.gr <- function(theta, args){
  n.theta <- length(theta)
  beta <- theta[1:(n.theta-1)]
  lnsigma <- theta[n.theta]
  
  y <- args$y
  X <- args$X
  Xb <- as.vector(X %*% beta)
  sigma <- exp(lnsigma)
  
  beta.gr <- apply((y-Xb)*X, 2, sum) / sigma**2
  lnsigma.gr <- sum((y-Xb)**2/sigma**2 - 1)
  return(-c(beta.gr, lnsigma.gr))
}
