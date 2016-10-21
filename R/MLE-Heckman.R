#' Heckman Models Estimation with MLE, Heckit
#' 
#' @name mle_heckman
#' @param y explained variable
#' @param X explanatory variable in matrix
#' @param z selection explained variable, binary
#' @param W selection explanatory variable in matrix
#' 
#' @return a list with 5 elements, mle result
#'   \item{method}{'mle'}
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table}
NULL

heckman.lnlike <- function(theta, args){
  # y~X * beta is the regression model
  # z~W * gamma is the selection model
  y <- args$y
  X <- args$X
  z <- args$z
  W <- args$W
  
  # theta = [gamma, theta, atanh(rho), ln(sigma)]
  n.beta <- dim(X)[2]
  n.gamma <- dim(W)[2]
  n.theta <- length(theta)
  gamma <- theta[1:n.gamma]
  beta <- theta[(n.gamma+1):(n.gamma+n.beta)]
  atanhrho <- theta[n.theta-1]
  lnsigma <- theta[n.theta]
  
  rho <- (exp(2*atanhrho)-1) / (exp(2*atanhrho)+1)
  sigma <- exp(lnsigma)
  
  Xb <- as.vector(X %*% beta)
  Wg <- as.vector(W %*% gamma)
  lnlike.z1 <- pnorm((Wg+(y-Xb)*rho/sigma)/sqrt(1-rho**2), log.p=TRUE) - 1/2*((y-Xb)/sigma)**2 - log(sqrt(2*pi)*sigma)
  lnlike.z0 <- pnorm(-Wg, log.p=TRUE)
  # Since 0*NA=NA, so do not use z*lnlike.z1 + (1-z)*lnlike.z0
  lnlike <- ifelse(z==1, lnlike.z1, lnlike.z0)
  return(-sum(lnlike))
}

heckman.gr <- function(theta, args){
  # y~X * beta is the regression model
  # z~W * gamma is the selection model
  y <- args$y
  X <- args$X
  z <- args$z
  W <- args$W
  
  # theta = [gamma, theta, atanh(rho), ln(sigma)]
  n.beta <- dim(X)[2]
  n.gamma <- dim(W)[2]
  n.theta <- length(theta)
  gamma <- theta[1:n.gamma]
  beta <- theta[(n.gamma+1):(n.gamma+n.beta)]
  atanhrho <- theta[n.theta-1]
  lnsigma <- theta[n.theta]
  
  rho <- (exp(2*atanhrho)-1) / (exp(2*atanhrho)+1)
  sigma <- exp(lnsigma)
  
  Xb <- as.vector(X %*% beta)
  Wg <- as.vector(W %*% gamma)
  nv <- (Wg+(y-Xb)*rho/sigma) / sqrt(1-rho**2)
  gamma.gr.z1 <- dnorm(nv)/pnorm(nv)/sqrt(1-rho**2) * W
  beta.gr.z1 <- (dnorm(nv)/pnorm(nv)*rho/sigma/sqrt(1-rho**2) - (y-Xb)/sigma**2) * (-X)
  atanhrho.gr.z1 <- dnorm(nv)/pnorm(nv)*((y-Xb)/sigma + rho*Wg)/(1-rho**2)**(3/2) * 4*exp(2*atanhrho)/(exp(2*atanhrho)+1)**2
  lnsigma.gr.z1 <- (dnorm(nv)/pnorm(nv)*(y-Xb)*rho/sqrt(1-rho**2)*(-1/sigma**2) + (y-Xb)**2/sigma**3 - 1/sigma) * sigma
  gamma.gr.z0 <- dnorm(-Wg)/pnorm(-Wg) * (-W)
  
  # ifelse(condition, r1, r2) return the result in shape of condition
  gamma.gr <- apply(ifelse(matrix(rep(z==1, n.gamma), ncol=n.gamma), gamma.gr.z1, gamma.gr.z0), 2, sum)
  beta.gr <- apply(ifelse(matrix(rep(z==1, n.beta), ncol=n.beta), beta.gr.z1, 0), 2, sum)
  lnsigma.gr <- sum(ifelse(z==1, lnsigma.gr.z1, 0))
  atanhrho.gr <- sum(ifelse(z==1, atanhrho.gr.z1, 0))
  return(-c(gamma.gr, beta.gr, atanhrho.gr, lnsigma.gr))
}

#' @export 
#' @rdname mle_heckman
mle.heckman.estimate <- function(y, X, z, W){
  args <- list(y=y, X=as.matrix(X), z=z, W=as.matrix(W))
  args <- clean4regression(args, assign.labels=c('z', 'W'), add.const=FALSE)
  args <- clean4regression(args, assign.labels=c('y', 'X'), exempt.cond=(args$z==0))
  
  params0 <- rep(1e-5, dim(args$W)[2]+dim(args$X)[2]+2)
  names(params0) <- c(paste0('W.', colnames(args$W)), 
                      paste0('X.', colnames(args$X)), 
                      'atanhrho', 'lnsigma')
  model.res <- mle.model(heckman.lnlike, args, params0=params0, gr=heckman.gr)
  return(model.res)
}
