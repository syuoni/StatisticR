#' Within Estimate
#' 
#' @name ols_within_estimate
#' @param y explained variable
#' @param X explanatory variable in matrix
#' @param group.indic group indicator variable
#' @param robust whether use robust convariance
#' 
#' @return a list with 8 elements
#'   \item{method}{'ols'}
#'   \item{convariance}{1 if robust, 0 if not}
#'   \item{observations}{number of samples}
#'   \item{R.square}{R^2}
#'   \item{adj.R.square}{adjusted R^2}
#'   \item{F.statistic}{F-statistic}
#'   \item{Prob.F}{Prob(F)}
#'   \item{table}{regression result table}
#' 
#' @export
ols.within.estimate <- function(y, X, group.indic, robust=FALSE){
  args <- list(y=y, X=as.matrix(X), group.indic=group.indic)
  args <- clean4regression(args, add.const=FALSE)
  y <- args$y
  X <- args$X
  group.indic <- args$group.indic
  
  for(lv in unique(group.indic)){
    y[group.indic==lv] <- y[group.indic==lv] - mean(y[group.indic==lv])
    # Propagation in R: two arrays are all raveled, 
    # then repeat the shorter one to meet the longer one
    if(dim(X)[2] > 1){
      X[group.indic==lv, ] <- t(t(X[group.indic==lv, ]) - apply(X[group.indic==lv, ], 2, mean))
    }else{
      X[group.indic==lv, ] <- X[group.indic==lv, ] - mean(X[group.indic==lv, ])
    }
  }
  model.res <- ols.estimate(y, X, robust=robust)
  return(model.res)
}
