#' Export ONE regression model result as a vector.
#' 
#' @name model_res_export
#' @param model.res        ONE regression model result
#' @param all.params       all potential parameters 
#' @param parentheses.type 't-statistic' or 'standard error' in parentheses
#' @param digits           digits
#' 
#' @return a vector of regression model result
#' 
#' @export 

model.res.export <- function(model.res, all.params=NULL, parentheses.type='t.statistic', digits=3){
  # if NO all.params specified, use all params from model.res
  if(is.null(all.params)){
    all.params <- row.names(model.res$table)
  }
  index <- NULL
  index[seq(1, by=2, along.with=all.params)] <- all.params
  index[seq(2, by=2, along.with=all.params)] <- paste0(all.params, '.par')
  
  row.index <- 1:dim(model.res$table)[1]
  names(row.index) <- row.names(model.res$table)
  row.index <- row.index[all.params]
  # row.index <- row.index[!is.na(row.index)]
  
  coef.str <- sapply(model.res$table[row.index, 'coef'], get.formatted, digits=digits)
  sig.stars <- sapply(model.res$table[row.index, 'p.value'], get.sig.stars)
  coef.stars <- paste0(coef.str, sig.stars)
  
  if(parentheses.type=='t.statistic'){
    in.parentheses <- sapply(model.res$table[row.index, 't.statistic'], get.formatted, digits=digits, parentheses=TRUE)
  }else{
    in.parentheses <- sapply(model.res$table[row.index, 'est.std.err'], get.formatted, digits=digits, parentheses=TRUE)
  }
  
  export.col <- NULL
  export.col[seq(1, by=2, along.with=all.params)] <- coef.stars
  export.col[seq(2, by=2, along.with=all.params)] <- in.parentheses
  
  # Additional Infomation
  if(model.res$method=='ols'){
    export.col <- c(export.col, model.res$observations, '', get.formatted(model.res$adj.R.square, digits=digits),
                    get.formatted(model.res$F.statistic, digits=digits))
  }else if(model.res$method=='mle'){
    export.col <- c(export.col, model.res$observations, get.formatted(model.res$lnlikelihood, digits=digits), '', '')
  }
  names(export.col) <- c(index, 'observations', 'lnlike', 'adj.R-square', 'F-statistic')
  return(export.col)
}

get.sig.stars <- function(p){
  if(is.na(p)){
    return('')
  }else if(p<=0.01){
    return('***')
  }else if(p<=0.05){
    return('**')
  }else if(p<=0.1){
    return('*')
  }else{
    return('')
  }
}

get.formatted <- function(v, digits=3, parentheses=FALSE){
  if(is.na(v)){
    return('')
  }
  v.str <- format(round(v, digits=digits), nsmall=digits, trim=TRUE)
  if(parentheses){
    return(paste0('(', v.str, ')'))
  }else{
    return(v.str)
  }
}
