#' Permute vector or matrix using a PBSmat object
#'
#' @description Permute vector or matrix using a PBSmat object
#'
#' @param x vector or matrix
#' @param PBSmat PBSmat object
#' @export
#'
PBS_perm <- function(x, PBSmat){UseMethod("PBS_perm")}

#' @export
PBS_perm.numeric <- function(x, PBSmat){
  type <- attr(PBSmat,"type")
  if(grepl("S",type)){
    sign <- sign (PBSmat)
  }else {sign<- 1L}
  PBSmat <- abs(PBSmat)
  PBSmat <- apply(PBSmat,2,function(pi){
    as.numeric(x[pi])
  })*sign
  unclass(PBSmat)
}

#' @export
PBS_perm.matrix <- function(x, PBSmat){
  type <- attr(PBSmat,"type")

  if(grepl("S",type)){
    PBSmat <- lapply(1:ncol(PBSmat),function(coli){x[abs(PBSmat[,coli]),]*sign(PBSmat[,coli]) })
  }else {
    PBSmat <- lapply(1:ncol(PBSmat),function(coli){x[PBSmat[,coli],] })

  }

  PBSmat
}
