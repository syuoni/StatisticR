#' @title Duration Models With Selection Estimation with MLE
#' @author syuoni
#' 
#' @references Boehmke F J, Morey D S, Shannon M. Selection Bias and Continuous‐Time Duration Models: 
#'   Consequences and a Proposed Solution[J]. American Journal of Political Science, 2006, 50(1): 192-207.
#' 
#' @name mle_duration_models_with_selection
#' @description Duration Models With Selection Estimation with MLE
#' @param t explained variable, vector, duration time
#' @param X explanatory variable in matrix
#' @param d explanatory variable, vector, whether the sample is failed (1=failure)
#' @param z selection explained variable, binary
#' @param W selection explanatory variable in matrix
#' 
#' @return a list with 4 elements, mle result
#'   \item{convergence}{1 if converge successfully, 0 if fail}
#'   \item{observations}{number of samples}
#'   \item{lnlikelihood}{log(likelihood)}
#'   \item{table}{regression result table}
NULL

exp.selection.lnlike <- function(theta, args){
  # t~X * beta is the duration model
  # z~W * gamma is the selection model
  # d=0 means right censored
  t <- args$t
  X <- args$X
  z <- args$z
  W <- args$W
  d <- args$d
  
  # theta = [gamma, theta, atanh(alpha)]
  n.beta <- dim(X)[2]
  n.gamma <- dim(W)[2]
  n.theta <- length(theta)
  gamma <- theta[1:n.gamma]
  beta <- theta[(n.gamma+1):(n.gamma+n.beta)]
  atanhalpha <- theta[n.theta]
  
  alpha <- (exp(2*atanhalpha)-1) / (exp(2*atanhalpha)+1)
  
  Xb <- as.vector(X %*% beta)
  Wg <- as.vector(W %*% gamma)
  lambda1 <- exp(-Wg)
  lambda2 <- exp(Xb)
  
  lnlike.z1d0 <- -lambda1-lambda2*t+log(lambda2)+log(1+alpha*(1-2*exp(-lambda2*t))*(1-exp(-lambda1)))
  lnlike.z1d1 <- -lambda1-lambda2*t+log(1+alpha*(1-exp(-lambda2*t))*(1-exp(-lambda1)))
  lnlike.z0 <- log(1-exp(-lambda1))
  lnlike <- ifelse(z==0, lnlike.z0, ifelse(d==0, lnlike.z1d0, lnlike.z1d1))
  return(-sum(lnlike))
}

#' @export 
#' @rdname mle_duration_models_with_selection
mle.exp.selection.estimate <- function(t, X, z, W, d){
  args <- list(t=t, X=as.matrix(X), z=z, W=as.matrix(W), d=d)
  args <- clean4regression(args, assign.labels=c('z', 'W'), add.const=FALSE)
  args <- clean4regression(args, assign.labels=c('t', 'X', 'd'), exempt.cond=(args$z==0))
  
  params0 <- rep(1e-5, dim(args$W)[2]+dim(args$X)[2]+1)
  names(params0) <- c(paste0('W.', colnames(args$W)),
                      paste0('X.', colnames(args$X)),
                      'atanhalpha')
  model.res <- mle.model(exp.selection.lnlike, args, params0=params0, gr=NULL)
  return(model.res)
}

weibull.selection.lnlike <- function(theta, args){
  # t~X * beta is the duration model
  # z~W * gamma is the selection model
  # d=0 means right censored
  t <- args$t
  X <- args$X
  z <- args$z
  W <- args$W
  d <- args$d
  
  # theta = [gamma, theta, atanh(alpha), ln(p)]
  n.beta <- dim(X)[2]
  n.gamma <- dim(W)[2]
  n.theta <- length(theta)
  gamma <- theta[1:n.gamma]
  beta <- theta[(n.gamma+1):(n.gamma+n.beta)]
  atanhalpha <- theta[n.theta-1]
  lnp <- theta[n.theta]
  
  alpha <- (exp(2*atanhalpha)-1) / (exp(2*atanhalpha)+1)
  p <- exp(lnp)
  
  Xb <- as.vector(X %*% beta)
  Wg <- as.vector(W %*% gamma)
  lambda1 <- exp(-Wg)
  lambda2 <- exp(Xb)

  lnlike.z1d0 <- -lambda1+log(1+alpha*(1-2*exp(-(lambda2*t)**p))*(1-exp(-lambda1)))+lnp+log(lambda2)+(p-1)*log(lambda2*t)-(lambda2*t)**p
  lnlike.z1d1 <- -lambda1-(lambda2*t)**p+log(1+alpha*(1-exp(-(lambda2*t)**p))*(1-exp(-lambda1)))
  # Stata/dursel/wblsel: Calculate lnlike.z1d1 with following formula
  # These two formulas are equivalent in theory, but a little diffrent numerically
  # lnlike.z1d1 <- log(1-(1-exp(-lambda1))-(1-exp(-(lambda2*t)**p))+(1-exp(-lambda1))*(1-exp(-(lambda2*t)**p))*(1+alpha*exp(-(lambda2*t)**p-lambda1)))
  lnlike.z0 <- log(1-exp(-lambda1))
  lnlike <- ifelse(z==0, lnlike.z0, ifelse(d==0, lnlike.z1d0, lnlike.z1d1))
  return(-sum(lnlike))
}

#' @export 
#' @rdname mle_duration_models_with_selection
mle.weibull.selection.estimate <- function(t, X, z, W, d){
  args <- list(t=t, X=as.matrix(X), z=z, W=as.matrix(W), d=d)
  args <- clean4regression(args, assign.labels=c('z', 'W'), add.const=FALSE)
  args <- clean4regression(args, assign.labels=c('t', 'X', 'd'), exempt.cond=(args$z==0))
  
  # params0 <- c(-.08711482, -.07519727, .01080475, .03833166, .02785502, .18414317, -.05496903, -1.7385206,
  #              .80768861, -1.9593325, .51412655, -.14661611, .3337918, .37420079, 1.0725059, -.65164642)
  # print(weibull.selection.lnlike(params0, args))
  params0 <- rep(1e-5, dim(args$W)[2]+dim(args$X)[2]+2)
  names(params0) <- c(paste0('W.', colnames(args$W)),
                      paste0('X.', colnames(args$X)),
                      'atanhalpha', 'lnp')
  model.res <- mle.model(weibull.selection.lnlike, args, params0=params0, gr=NULL)
  return(model.res)
}
