#' create permutation, bootstrap or shuffle matrix
#'
#' @description create permutation, bootstrap or shuffle matrix
#'
#' @param n integer, length of vector to transform
#' @param np integer, number of tranformations
#' @param type character intecating the type of transormation. "P", "S", "B", "PS", "BS" are available.
#' @importFrom stats rbinom
#' @export
PBSmat = function(n, np = 4000, type = "P"){
  if(type == "SP"){type="PS"}
  if(type == "SB"){type="BS"}

  ##
  if(grepl(pattern = "P",type)){
    mat <- replicate(np - 1, sample(n, n, replace = F))}
  else if(grepl(pattern = "B",type)){
    mat <- replicate(np - 1, sample(n, n, replace = T))}
  else if(type=="S"){
    mat <- replicate(np-1, 1:n)}

  if(grepl(pattern = "S",type)){
    neg = rbinom(n*(np-1),1,prob=0.5)==1
    mat[neg] = -(mat[neg])
  }

  mat = cbind(1:n,mat)
  attr(mat,"type") = type
  class(mat)=c("PBSmat",class(mat))
  mat
}
