#' Compute the unmodified blup of a lmerMod object
#'
#' @description Compute the unmodified blup of a lmerMod object
#'
#' @param model a lmerMod object.
#' @export
blup_blup <- function(model){
  Uhat.list <- ranef(model)
  Uhat.list <- lapply(Uhat.list,as.matrix.data.frame)
  ehat <- resid(model)
  out = list(Uhat_list = Uhat.list, ehat = ehat)
  return(out)

}
