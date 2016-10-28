#' Generate dummy variables from ONE variable
#' 
#' @name gen_dummy
#' @param vec    variable vector
#' @param levels variable values, each generates a dummy
#' @param prefix prefix for dummy variables
#' @param head   dummy variable names
#' @param bool   return bool dummy variables if TRUE
#' 
#' @return dummy variable matrix
#' 
#' @export 
gen.dummy <- function(vec, levels=NULL, prefix=NULL, head=NULL, bool=FALSE){
  if(is.null(levels)){
    levels <- unique(vec)
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
  return(dummies)
}
