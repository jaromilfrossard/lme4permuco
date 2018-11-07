#' Compute the contrasts of the unmodified blup of a lmerMod object
#'
#' @description Compute the contrasts unmodified blup of a lmerMod object
#'
#' @param model a lmerMod object.
#' @export
blup_contrblup <- function(model){
  Uhat.list <- ranef(model)
  Uhat.list <- lapply(Uhat.list,as.matrix.data.frame)
  contrlist <- getContrlist(model)

  #length of sampling unit for full contrasts
  SUn <- lapply(names(model@reTrms$contrlist),function(namei) {
    namei = unlist(strsplit(unlist(strsplit(namei, "[|]"))[2],"[:]"))[1]
    length(levels(droplevels(model@frame[,gsub(" ", "", namei, fixed = TRUE)])))

  })



  Uhat.list <- mapply(function(ni,contri,uhati)(Diagonal(ni)%x%contri)%*%uhati,
                      contri=contrlist,ni=SUn, uhati = Uhat.list,SIMPLIFY = F)




  ehat <- resid(model)
  out = list(Uhat_list = Uhat.list, ehat = ehat)
  return(out)

}
