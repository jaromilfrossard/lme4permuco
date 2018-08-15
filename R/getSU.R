##get sampling unit after find bars
getSU <- function(bars){
  if(length(bars)==1){
    return(bars)
  }
  else{
    getSU(bars[[2]])
  }
}
