#' @title Clean for Regression
#' @author syuoni
#' 
#' @name clean4regression
#' @description drop NA observations and add constant variable, prepare for regression model
#' @param args          arguments for regression model
#' @param assign.labels NA samples in args with specific labels would be dropped 
#' @param exempt.cond   vector, exempt condition, sample would not be dropped if it is TRUE
#' @param add.const     whether add const variable for matrix
#' 
#' @return a list, args for regression
#' 
#' @export 

clean4regression <- function(args, assign.labels=NULL, exempt.cond=FALSE, add.const=TRUE){
  if(is.null(assign.labels)){
    assign.labels <- names(args)
  }
  na.indicator <- NULL
  for(lb in assign.labels){
    na.indicator <- cbind(na.indicator, is.na(args[[lb]]))
  }
  drop.indicator <- apply(na.indicator, 1, any) & (!exempt.cond)
  
  model.args <- list()
  for(lb in names(args)){
    if(is.vector(args[[lb]])){
      model.args[[lb]] <- args[[lb]][!drop.indicator]
    }else{
      model.args[[lb]] <- args[[lb]][!drop.indicator, ]
      
      # for matrix: drop variables with unique value
      # so if const variable is in matrix, it would be dropped, but may be added latter if add.const=TRUE
      model.args[[lb]] <- model.args[[lb]][, apply(model.args[[lb]], 2, function(vec){
        if(length(unique(vec)) > 1){
          return(TRUE)
        }else{
          return(FALSE)
        }
      })]

      if(add.const){
        model.args[[lb]] <- cbind(model.args[[lb]], `_const`=1)
      }
    }
  }
  return(model.args)
}
