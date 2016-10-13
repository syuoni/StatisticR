#' @title model results export
#' @author syuoni
#' 
#' @name mle_model_res_export
NULL

#' @description export ONE mle model as a vector
#' @param model.res ONE mle regression model result
#' 
#' @return a vector
#' @export
#' @rdname mle_model_res_export

mle.res.export <- function(model.res, ...){
  export.col <- model.res.export(model.res, ...)
  index <- names(export.col)
  export.col <- c(export.col, model.res$observations, get.formatted(model.res$lnlikelihood))
  names(export.col) <- c(index, 'observations', 'lnlike')
  return(export.col)
}

#' @description export MULTIPLE mle models as a table (data frame)
#' @param model.res.list mle regression models results
#' 
#' @return a data frame
#' @export
#' @rdname mle_model_res_export
mle.res.table.export <- function(model.res.list, ...){
  res.table <- sapply(model.res.list, mle.res.export, ...)
  res.table <- data.frame(res.table)
  colnames(res.table) <- paste0('(', 1:length(res.table), ')')
  return(res.table)
}