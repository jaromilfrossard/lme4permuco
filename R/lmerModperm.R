#' Permutation test on lmerMod object
#'
#' @description permutation test on lmerMod object
#'
#' @param model a lmerMod object
#' @param blup_FUN the type of BLUP to permute. Default is CGR.
#' @param ... other arguments
#' @export
lmerModperm <- function(model, blup_FUN, ...){UseMethod("lmerModperm")}

#' @export
lmerModperm.lmerModgANOVA <- function(model, blup_FUN = blup_cgr, ...){

  blup <- blup_FUN(model)

  out <- list()
  out$model <- model
  out$blup <- blup
  out

}
