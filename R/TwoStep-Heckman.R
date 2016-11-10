#' Heckman Models Estimation with Two Step, i.e. Heckit
#' 
#' @name twostep_heckman
#' @param y explained variable
#' @param X explanatory variable in matrix
#' @param z selection explained variable, binary
#' @param W selection explanatory variable in matrix
#' 
#' @return a list with 5 elements
#'   \item{method}{'two-step'}
#'   \item{selection}{selection model result of probit}
#'   \item{regression}{regression model result of ols}
#'   \item{sigma}{sigma}
#'   \item{rho}{rho}
#' 
#' @export 
twostep.heckman.estimate <- function(y, X, z, W){
  args <- list(y=y, X=as.matrix(X), z=z, W=as.matrix(W))
  args <- clean4regression(args, assign.labels=c('z', 'W'), add.const=FALSE)
  args <- clean4regression(args, assign.labels=c('y', 'X'), exempt.cond=(args$z==0), add.const=FALSE)
  y <- args$y
  X <- args$X
  z <- args$z
  W <- args$W
  
  # Selection model
  selection.res <- mle.probit.estimate(z, W)
  
  # Inverse Mills Ratio, which is denoted as lambda
  z.star <- cbind(W, 1) %*% selection.res$table$coef
  lambda <- dnorm(-z.star)/(1-pnorm(-z.star))
  lambda.z1 <- as.vector(lambda)[z==1]
  
  # Regression Model
  y.z1 <- y[z==1]
  # Note if X is reduced to vector, cbind would always transform the final result to matrix. 
  X.z1 <- cbind(X[z==1, ], lambda.z1)
  colnames(X.z1) <- c(colnames(X), 'lambda')
  regression.res <- ols.estimate(y.z1, X.z1, robust=TRUE)
  
  # Calculation for rho and sigma, where rho is the correlation for the 2 error terms,
  # and sigma is the variance for regression model error term. 
  # This calculation result is verified by an example from Stata, but the formula remains to be confirm. 
  rho.sigma.hat <- regression.res$table['lambda', 'coef']
  residual.z1 <- y.z1 - cbind(X.z1, 1) %*% regression.res$table$coef
  sigma.hat <- sqrt(mean(residual.z1**2) + mean(lambda.z1*(z.star[z==1]+lambda.z1)) * rho.sigma.hat**2)
  rho.hat <- rho.sigma.hat / sigma.hat
  
  model.res <- list(method      ='two-step',
                    selection   =selection.res,
                    regression  =regression.res,
                    sigma       =sigma.hat,
                    rho         =rho.hat)
  return(model.res)
}
