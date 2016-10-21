#' export MULTIPLE models as a table (data frame)
#' 
#' @name model_res_table_export
#' @param model.res.list a list of several regression model results
#' @param all.params     all potential parameters 
#' @param ...            parameters passed to function model_res_export
#' 
#' @return a data frame
#' @export
model.res.table.export <- function(model.res.list, all.params=NULL, ...){
  # if NO all.params specified, use all params from model.res.list
  if(is.null(all.params)){
    all.params <- NULL
    for(model.res in model.res.list){
      all.params <- c(all.params, row.names(model.res$table))
    }
    all.params <- unique(all.params)
  }
  
  res.table <- sapply(model.res.list, model.res.export, all.params=all.params, ...)
  res.table <- data.frame(res.table)
  colnames(res.table) <- paste0('(', 1:length(res.table), ')')
  return(res.table)
}