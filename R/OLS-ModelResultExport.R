#' @title model results export for OLS
#' @author syuoni
#' 
#' @name ols_model_res_export
NULL

#' @description export ONE ols model as a vector
#' @param model.res ONE ols regression model result
#' @param ... parameters passed
#' 
#' @return a vector
#' @export
#' @rdname ols_model_res_export
ols.res.export <- function(model.res, ...){
  export.col <- model.res.export(model.res, ...)
  index <- names(export.col)
  export.col <- c(export.col, model.res$observations, get.formatted(model.res$adj.R.square),
                  get.formatted(model.res$F.statistic))
  names(export.col) <- c(index, 'observations', 'R-square', 'F-statistic')
  return(export.col)
}

