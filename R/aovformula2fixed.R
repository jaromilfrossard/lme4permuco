# #'@param formula a aov type formula
# #'@param n_SU number of random unit
aovformula2fixed = function(formula, n_SU){
  if(n_SU==1){
    f = formula(paste(deparse(formula[[2]]),"~",deparse(formula[[3]][[2]]),sep=""))
    return(f)
  }
  else{
    formula[[3]] <- formula[[3]][[2]]
    return(aovformula2fixed(formula = formula,n_SU = n_SU - 1))
  }
}
