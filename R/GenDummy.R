#' Generate dummy variables from ONE variable
#' 
#' @name gen_dummy
#' @param vec    variable vector
#' @param levels variable values, each generates a dummy
#' @param head   dummy variable names
#' 
#' @return dummy variable matrix
#' 
#' @export 
gen.dummy <- function(vec, levels=NULL, head=NULL){
  if(is.null(levels)){
    levels <- unique(vec)
  }
  dummies <- NULL
  for(lv in levels){
    dummies <- cbind(dummies, as.integer(vec==lv))
  }
  colnames(dummies) <- head
  return(dummies)
}
