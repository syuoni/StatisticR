#' @title model results export
#' @author syuoni
#' 
#' @description Export ONE regression model result as a vector
#' 
#' @name model_res_export
#' @param model.res        ONE regression model result
#' @param all.params       all potential parameters 
#' @param parentheses.type 't-statistic' or 'standard error'
#' @param digits
#' 
#' @return a vector
#' 
#' @export 

model.res.export <- function(model.res, all.params, parentheses.type='t.statistic', digits=3){
  index <- NULL
  index[seq(1, by=2, along.with=all.params)] <- all.params
  index[seq(2, by=2, along.with=all.params)] <- paste0(all.params, '.par')
  
  coef.str <- sapply(model.res$table[all.params, 'coef'], get.formatted)
  sig.stars <- sapply(model.res$table[all.params, 'p.value'], get.sig.stars)
  coef.stars <- paste0(coef.str, sig.stars)
  
  if(parentheses.type=='t.statistic'){
    in.parentheses <- sapply(model.res$table[all.params, 't.statistic'], get.formatted, parentheses=TRUE)
  }else{
    in.parentheses <- sapply(model.res$table[all.params, 'est.std.err'], get.formatted, parentheses=TRUE)
  }
  
  export.col <- NULL
  export.col[seq(1, by=2, along.with=all.params)] <- coef.stars
  export.col[seq(2, by=2, along.with=all.params)] <- in.parentheses
  names(export.col) <- index
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
