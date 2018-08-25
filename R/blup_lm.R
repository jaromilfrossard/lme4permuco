#' Compute the blup as fixed effect of a lmerMod object
#'
#' @description Compute the blup as fixed effect of a lmerMod object
#'
#' @param model a lmerMod object.
#' @export
blup_lm <- function(model){
  UseMethod("blup_lm")
}


blup_lm.lmerModgANOVA = function(model){

  X = getX(model)

  Ztlist = getZtlm(model,"reduced")
  contrlist = getZtlm(model,"contrasts")

  Zt = do.call("rbind",Ztlist)
  qr_xz=qr(cbind(Matrix(X,sparse=T),t(Zt)))

  beta = qr.coef(qr_xz,getME(model,"y"))

  gamma = beta[-c(1:ncol(X ))]
  assign = unlist(sapply(1:length(Ztlist), function(i)rep(i,nrow(Ztlist[[i]]))))


  gamma = lapply(1:length(Ztlist), function(i)matrix(gamma[assign==i],ncol=1))
  names(gamma) = names(Ztlist)

  gamma = mapply(function(gi,ci){ci%*%gi},
         gi = gamma, ci = contrlist,
         SIMPLIFY = F
         )
  ehat = qr.resid(qr_xz,getME(model,"y"))

  out = list(Uhat_list = gamma, ehat = matrix(ehat,ncol=1))
  return(out)



}



blup_lm.list = function(model){
  if(sum(names(model)== c("fr","X","reTrms","REML","formula","wmsgs" ))!=6){
    stop("the model is not the output of gANOVA_lFormula()")
  }

  X = getX(model)
  formula = as.formula(model$formula)
  formula_aov = lme4formula2aov(formula)

  Ztlist = getZtlm(formula_aov,"reduced",data=model$fr)
  contrlist = getZtlm(formula_aov,"contrasts",data=model$fr)

  Zt = do.call("rbind",Ztlist)
  qr_xz=qr(cbind(Matrix(X,sparse=T),t(Zt)))

  beta = qr.coef(qr_xz,model$fr[,attr(terms(model$fr),"response")])

  gamma = beta[-c(1:ncol(X ))]
  assign = unlist(sapply(1:length(Ztlist), function(i)rep(i,nrow(Ztlist[[i]]))))


  gamma = lapply(1:length(Ztlist), function(i)matrix(gamma[assign==i],ncol=1))
  names(gamma) = names(Ztlist)

  gamma = mapply(function(gi,ci){ci%*%gi},
                 gi = gamma, ci = contrlist,
                 SIMPLIFY = F
  )
  ehat = qr.resid(qr_xz,model$fr[,attr(terms(model$fr),"response")])

  out = list(Uhat_list = gamma, ehat = matrix(ehat,ncol=1))
  return(out)



}
