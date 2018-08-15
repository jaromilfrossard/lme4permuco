#' Extract contrast of the random effect for a gANOVA model
#'
#' @description Extract the contrast (polynomial) of the random effect and with their associated matrix and coding (R contr.poly)
#'
#' @param model a gANOVA model
#' @importFrom stats contr.poly formula model.matrix terms formula
#' @importFrom lme4 findbars factorize
#' @export
ranefcontrasts <- function(model){UseMethod("ranefcontrasts")}

#' @export
ranefcontrasts.lmerModgANOVA <- function(model){
  formula <- formula(model)
  bars <- findbars(lme4:::RHSForm(terms(formula)))
  names(bars) <- lme4:::barnames(bars)
  fr <-  model@frame
  rf <- ranef(model)

  bars <- bars[match(names(bars),names(rf))]

  contrlist <- lapply(bars, mkcontrastslist, fr)

  ranefcontr = list()
  for(i in 1:length(contrlist)){
    ranefcontr[[i]]<- as.numeric(contrlist[[i]]$contr%*%matrix(as.numeric(as.matrix(rf[[i]])),ncol=contrlist[[i]]$n))}
  names(ranefcontr) <- names(contrlist)

  out = list()
  out$ranef = ranefcontr
  out$contr = lapply(contrlist,function(x)x$contr)
  out$Zmat = lapply(contrlist,function(x)x$sm)
  out

}
