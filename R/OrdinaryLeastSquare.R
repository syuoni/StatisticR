#' @title Ordinary Least Square
#' @author syuoni
#' 
#' @name ordinary_least_square
#' @description Ordinary Least Square
#' @param y explained variable
#' @param X explanatory variable in matrix
#' @param robust whether use robust convariance
#' 
#' @return a list with 7 elements
#'   \item{convariance}{1 if robust, 0 if not}
#'   \item{observations}{number of samples}
#'   \item{R.square}{R^2}
#'   \item{adj.R.square}{adjusted R^2}
#'   \item{F.statistic}{F-statistic}
#'   \item{Prob.F}{Prob(F)}
#'   \item{table}{regression result table}
#' 
#' @export 

ols.estimate <- function(y, X, robust=FALSE){
  args <- list(y=y, X=as.matrix(X))
  args <- clean4regression(args)
  y <- args$y
  X <- args$X
  n <- dim(X)[1]
  K <- dim(X)[2]
  
  XtX <- t(X) %*% X
  XtXi <- solve(XtX)
  # b.hat = B %*% y
  B <- XtXi %*% t(X)
  # projection matrix
  P <- X %*% B
  # annihilator matrix
  M <- diag(n) - P
  
  b.hat <- B %*% y
  y.hat <- P %*% y
  e.hat <- M %*% y
  s.sqr <- sum(e.hat**2) / (n-K)
  s <- s.sqr ** 0.5
  
  # R-square
  y.mean <- mean(y)
  SST <- sum((y-y.mean)**2)
  SSE <- sum((y.hat-y.mean)**2)
  SSR <- sum(e.hat**2)
  R.sqr <- SSE / SST
  adj.R.sqr <- 1 - (SSR/(n-K)) / (SST/(n-1))
  
  if(robust==FALSE){
    convariance <- 'non-robust'
    # stanard errors and t-test
    est.cov <- s.sqr * XtXi
    est.std.err <- diag(est.cov) ** 0.5
    t.statistic <- b.hat / est.std.err
    p.value <- (1 - pt(abs(t.statistic), n-K)) * 2
    
    t95 <- qt(0.975, n-K)
    conf.int.lower <- b.hat - t95 * est.std.err
    conf.int.upper <- b.hat + t95 * est.std.err
    
    # F-test for all beta=0 (except for the constant)
    R <- diag(K)[1:K-1, ]
    Rb <- R %*% b.hat
    RXtXiRti <- solve(R %*% XtXi %*% t(R))
    F.statistic <- t(Rb) %*% RXtXiRti %*% Rb / ((K-1)*s.sqr)
    # matrix %*% matrix, or matrix %*% vector, or vector %*% vector, would always return matrix
    # while matrix[1, ] or matrix[, 1], it would be reduced to vector automatically
    F.statistic <- F.statistic[1]
    F.p.value <- 1 - pf(F.statistic, K-1, n-K)
  }else{
    convariance <- 'robust'
    # sqrt(n)*(b-beta) ~ N(0, Avarb)
    # b-beta = (X'X)-1 * X' * epsilon
    Avarb.hat <- n * B %*% diag(as.vector(e.hat)**2) %*% t(B)
    # freedom adjusted as n-K to be consistent with stata
    Avarb.hat <- Avarb.hat * n / (n-K)
    est.cov <- Avarb.hat / n
    est.std.err <- diag(est.cov) ** 0.5
    t.statistic <- b.hat / est.std.err
    
    # use t(n-K) distribution to be consistent with stata
    # option: use standard normal distribution
    p.value <- (1 - pt(abs(t.statistic), n-K)) * 2
    
    t95 <- qt(0.975, n-K)
    conf.int.lower <- b.hat - t95 * est.std.err
    conf.int.upper <- b.hat + t95 * est.std.err
    
    # Wald-test and F-test are equivalent in large sample
    # F-test for all beta=0 (except for the constant)
    R <- diag(K)[1:K-1, ]
    Rb <- R %*% b.hat
    
    RXtXiRt <- solve(R %*% XtXi %*% t(R))
    F.statistic <- t(Rb) %*% RXtXiRt %*% Rb / ((K-1)*s.sqr)
    
    RAvarb.hatRti <- solve(R %*% Avarb.hat %*% t(R))
    Wald.statistic <- n * t(Rb) %*% RAvarb.hatRti %*% Rb
    F.statistic <- Wald.statistic / (K-1)
    F.statistic <- F.statistic[1]
    F.p.value <- 1 - pf(F.statistic, K-1, n-K)
  }
  
  coef.table <- data.frame(coef=b.hat, 
                           est.std.err, t.statistic, p.value, 
                           conf.int.lower, conf.int.upper,
                           row.names=colnames(X))
  model.res <- list(convariance =convariance,
                    observations=n,
                    R.square    =R.sqr,
                    adj.R.square=adj.R.sqr,
                    F.statistic =F.statistic,
                    Prob.F      =F.p.value,
                    table       =coef.table)
  return(model.res)
}



