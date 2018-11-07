#' Extract a list of contrasts of the random design from gANOVA model
#'
#'@description Extract a list of contrasts of the random design from gANOVA model
#'@param model a lmerModgANOVA model

#'@export

getContrlist <- function(model){

  contrlist <- model@reTrms$contrlist
  contrlist <- lapply(contrlist,function(contrlisti){
    if(is.null(contrlisti)){
      contrlisti = Matrix(1,ncol = 1,nrow=1)
    }else{
      if(length(contrlisti)==1){
        contrlisti = contrlisti[[1]]}
      else{contrlisti = do.call("%x%",contrlisti)}}
    contrlisti
  })

  return(contrlist)

}
