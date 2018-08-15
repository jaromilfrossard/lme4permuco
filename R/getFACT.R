#get factors after find bars

getFACT <- function(bars){
  if(length(bars)==1){
    return(NULL)
  }
  else{
    c(bars[[3]], getFACT(bars[[2]]))
  }
}
