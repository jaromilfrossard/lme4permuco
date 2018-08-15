#' Transform a lme4 formula into an aov-like formula
#'
#' @description Transform a lme4 formula into a aov-like formula
#'
#' @param x a lmerMod or formula object.
#' @importFrom stats getCall
#' @export
lme4formula2aov <- function(x) UseMethod("lme4formula2aov")

#'@export
lme4formula2aov.lmerMod <- function(x){
  call = getCall(x)
  formula = eval(call$formula)
  lme4formula2aov.formula(formula)


}

#'@export
lme4formula2aov.formula <- function(x){
  fformula = lme4:::getFixedFormula(x)
  bars = findbars(x)
  SU = lapply(bars,function(b){
    getSU(b[[3]])})

  fact = lapply(bars,function(b){
    getFACT(b[[3]])})
  nullf = sapply(fact,is.null)
  #fact[nullf][[1]] = 1


  fact = lapply(fact,function(fl){
    paste(unlist(fl,recursive = T),collapse=":",sep=":")
  })
  uniqueSU = unique(SU)


  formula_aov=list()
  for(i in 1:length(uniqueSU)){
    sui = uniqueSU[[i]]
    selected_factor = unlist(fact[SU==sui])
    selected_factor=selected_factor[selected_factor!=""]

    if( length(selected_factor) == 0 ){
      formula_aov[[i]] = paste("+Error(",sui,")",collapse="")
    }else{

      formula_aov[[i]] = paste("+Error(",sui,"/(",
                               paste(unique(selected_factor),
                                     collapse ="+"),"))",sep="",collapse="")}
  }

  formula(paste(do.call("c",c(deparse(fformula),formula_aov)),collapse =""))
}
