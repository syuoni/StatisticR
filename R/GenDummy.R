#' Generate dummy variables from ONE variable
#' 
#' @name gen_dummy
#' @param vec      variable vector
#' @param levels   variable values, each generates a dummy
#' @param prefix   prefix for dummy variables
#' @param head     dummy variable names
#' @param first.rm remove the first dummy (as base) if TRUE
#' @param bool     return bool dummy variables if TRUE
#' 
#' @return dummy variable matrix
#' 
#' @export 
gen.dummy <- function(vec, levels=NULL, prefix=NULL, head=NULL, first.rm=FALSE, bool=FALSE){
  if(is.null(levels)){
    levels <- sort(unique(vec))
  }
  dummies <- NULL
  for(lv in levels){
    if(bool){
      dummies <- cbind(dummies, vec==lv)
    }else{
      dummies <- cbind(dummies, as.integer(vec==lv))
    }
  }
  if(!is.null(head)){
    colnames(dummies) <- head
  }else if(!is.null(prefix)){
    colnames(dummies) <- paste(prefix, levels, sep='.')
  }
  if(first.rm){
    if(dim(dummies)[2] > 2){
      dummies <- dummies[, 2:dim(dummies)[2]]
    }else{
      head <- colnames(dummies)
      dummies <- matrix(dummies[, 2])
      colnames(dummies) <- head[2]
    }
  }
  return(dummies)
}
