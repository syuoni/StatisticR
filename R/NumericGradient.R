#' @title Numeric Gradient Function
#' @author syuoni
#' 
#' @name numeric_gradient
#' @param func function to calculate numeric gradient function
#' 
#' @return gradient function
#' 
#' @export 

numeric.gr <- function(func){
  dx <- 1e-3
  func.gr <- function(theta, ...){
    n.theta <- length(theta)
    gr <- NULL
    for (i in 1:n.theta){
      theta0 <- theta
      theta1 <- theta
      theta0[i] <- theta[i] - dx/2
      theta1[i] <- theta[i] + dx/2
      gr <- c(gr, (func(theta1, ...)-func(theta0, ...)) / dx)
    }
    return(gr)
  }
  return(func.gr)
}

