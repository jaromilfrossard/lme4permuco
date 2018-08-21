#' Get the type of variable
#'
#' @description get the type of variable of a fixed design given a sampling unit
#'
#' @param formula a formula of a fixed design
#' @param data a dataframe
#' @SU a sampling unit
#' @export
getBetweenVar <- function(formula,data,SU){
  data = droplevels(data)
  termf <- terms(formula)
  termf <- delete.response(termf)
  ftable <- attr(termf,"factors")
  vars <-rownames(ftable)

  #return
  selvar = sapply(vars,function(vari){

    tabi = xtabs(formula(paste("~",vari,"+",SU,sep="",collapse="")),data=data)!=0
    if(sum(colSums(tabi)!=1)==0){
      T
    }else{
      F
    }})

  colSums(ftable [selvar,,drop=F])!=0

}
