#' Generate dummy variables from ONE variable
#' 
#' @name gen_dummy
#' @param vec     variable vector
#' @param levels  variable values, each generates a dummy
#' @param prefix  prefix for dummy variables
#' @param head    dummy variable names
#' @param last.rm remove the last dummy if TRUE
#' @param bool    return bool dummy variables if TRUE
#' 
#' @return dummy variable matrix
#' 
#' @export 
gen.dummy <- function(vec, levels=NULL, prefix=NULL, head=NULL, last.rm=FALSE, bool=FALSE){
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
  if(last.rm){
    if(dim(dummies)[2] > 2){
      dummies <- dummies[, 1:(dim(dummies)[2]-1)]
    }else{
      head <- colnames(dummies)
      dummies <- matrix(dummies[, 1])
      colnames(dummies) <- head[1]
    }
  }
  return(dummies)
}
