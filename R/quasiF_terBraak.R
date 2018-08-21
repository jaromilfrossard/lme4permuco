quasiF_terBraak <- function(args){
  data = args$model@frame
  formula <- formula(args$model)
  bars <- findbars(lme4:::RHSForm(formula))
  names(bars) <- lme4:::barnames(bars)
  SU = lapply(bars,function(b){
    getSU(b[[3]])})
  formula_fix <- lme4:::getFixedFormula(formula)
  SU = unique(sapply(SU, deparse))


  btable <- t(sapply(SU,function(sui){
    getBetweenVar(formula = formula_fix,data = data,SU = sui)
  }))

  formula_full = paste("~",paste(colnames(btable),collapse="+"),sep="")

  formula_full = paste(c(formula_full ,paste("Error(",sapply(1:nrow(btable),function(i){
    paste(rownames(btable)[i],"/(",paste(colnames(btable)[!btable[i,]],collapse="+"),")")
  }),")")),collapse=" + ")
  formula = formula(paste(deparse(formula[[2]]),formula_full))
  Ztlist = getZtlm(formula,"reduced",data=data)
}
