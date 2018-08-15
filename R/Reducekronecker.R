Reducekronecker <- function(list, ...) {
  if(length(list)==1){
    return(list[[1]])
  }else{
    list[[1]] = kronecker(list[[1]],list[[2]],...=...)
    list=list[-2]
    Reducekronecker(list=list,... = ...)
  }
  }
