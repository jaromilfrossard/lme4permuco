#' Get one or multiple permutation/bootstrap/shuffle from a pbsmat object
#'
#' @description Get one or multiple permutation/bootstrap/shuffle from a pbsmat object
#'
#' @param x a PBSmat object.
#' @param j the indicator of the transformation to get.
#' @export


getpbs <- function(x, j){ UseMethod("getpbs")}

getpbs.PBSmat <- function(x, j){
  type <- attr(x,"type")
  x <- matrix(x[ ,j, drop = F], nrow = nrow(x))
  attr(x,"type") <- type
  class(x) <- c("PBSmat",class(x))
  x
}
