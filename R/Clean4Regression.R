#' Clean data frame for fegression, including dropping NA observations 
#'   and adding constant variable, make it prepared for regression models
#' 
#' @name clean4regression
#' @param args          arguments for regression model, before cleaning
#' @param assign.labels NA samples in args with specific labels would be dropped 
#' @param exempt.cond   vector, exempt condition, sample would not be dropped if it is TRUE
#' @param add.const     whether add const variable for matrix
#' 
#' @return a list, arguments for regression model, after cleaning
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
  keep.nrow <- sum(!drop.indicator)
  
  model.args <- list()
  for(lb in names(args)){
    if(is.vector(args[[lb]])){
      model.args[[lb]] <- args[[lb]][!drop.indicator]
    }else{
      model.args[[lb]] <- args[[lb]][!drop.indicator, ]
      if(dim(args[[lb]])[2] == 1){
        # if the matrix has ONE column initially, we do NOT check whether this column is unique
        # Actually if it is unique, the model must be wrong
        model.args[[lb]] <- matrix(model.args[[lb]], ncol=1)
        colnames(model.args[[lb]]) <- colnames(args[[lb]])
      }else{
        # for matrix: drop variables with unique value
        # so if const variable is in matrix, it would be dropped, but may be added latter if add.const=TRUE
        not.unique <- apply(model.args[[lb]], 2, function(vec){
          if(length(unique(vec)) > 1){
            return(TRUE)
          }else{
            return(FALSE)
          }
        })
        # if just ONE column left, the matrix would be reduced to vector automatically
        # so we should transform it back to matrix manually
        model.args[[lb]] <- matrix(model.args[[lb]][, not.unique], nrow=keep.nrow)
        colnames(model.args[[lb]]) <- colnames(args[[lb]])[not.unique]
      }
      if(add.const){
        model.args[[lb]] <- cbind(model.args[[lb]], `_const`=1)
      }
    }
  }
  return(model.args)
}
