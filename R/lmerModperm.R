#' Permutation test on lmerMod object
#'
#' @description permutation test on lmerMod object
#'
#' @param model a lmerMod object
#' @param blupstar a character string indicating the type of random effect to permute. Default is "cgr" and "blup" and "lm"
#' @param np an integer indicating the number of permutation/bootstrap/suffle.
#' @param method a character string indicating the permutation/bootstrap/suffle method, default is "terBraak".
#' @param assigni a integer indicating the effect to test. Default is 1.
#' @param statistic a character string indicating the test statistic. default is "Satterthwaite".
#' @param ... other arguments
#' @importFrom stats pf anova as.formula
#' @export
lmerModperm <- function(model, blupstar, np, method, assigni, statistic, ...){UseMethod("lmerModperm")}

#' @export
lmerModperm.lmerModgANOVA <- function(model, blupstar = "cgr", np = 4000, method = "terBraak", assigni = 1, statistic = "Satterthwaite",...){
  argslist <- formals()
  mc <- match.call()
  argslist[names(as.list(mc[-1]))] = as.list(mc[-1])

  switch(blupstar,
        "cgr" = {blup_FUN = blup_cgr},
        "lm" = {blup_FUN = blup_lm},
        "blup" = {blup_FUN = blup_blup})


  switch(statistic,
         "Satterthwaite" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,5]}},
         "Satterthwaite_p" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,6]}},
         "Satterthwaite_logp" = {FUN_stat = function(model,assigni){
           ano = anova(model,type=3,ddf="Satterthwaite")
           abs(pf(q=ano[assigni,5], df1 = ano[assigni,3], df2 = ano[assigni,4], lower.tail = F, log.p = T))}},
         "quasiF" = {FUN_stat =function(model,assigni){model[,1]}},
         "quasiF_logp" = {FUN_stat =function(model,assigni){abs(pf(q = model[,1] , df1 = model[,10], df2 = model[,11],lower.tail = T, log.p = T))}}
         )


  if(statistic%in% c("quasiF","quasiF_logp") ){
    switch(method,
           "terBraak" = {FUN_p = function(x)quasiF_terBraak(x)},
           "dekker" = {stop("the dekker method do not work with quasi F statistics.")})
  }else{
    switch(method,
           "terBraak" = {FUN_p = function(x)lmerModperm_terBraak(x)},
           "dekker" = {FUN_p = function(x)lmerModperm_dekker(x)})
  }


  ##computing the blupstar
  blup <- blup_FUN(model)


  ### create new error

  reff <- blup$Uhat_list
  reff$error <- matrix(blup$ehat,ncol=1)



  PBSlist <- list(
    fixeff = PBSmat(n = length(getME(model,"y")), np = np, type = "P"),
    ranef = lapply(lapply(reff,length),function(n){
    PBSmat(n=n,np =np,type = "PS")
  }))


  if(eval(argslist$blupstar)=="lm"){
    Ztlist <- getZtlm(model,"reduced")
    X <- getME(model,"X")
    beta <- qr.coef(qr(cbind(X,t(do.call("rbind",Ztlist)))),getME(model,"y"))[1:ncol(X)]
    Ztlist <- getZtlm(model,"full")

  }else if( (eval(argslist$blupstar)=="blup")|
      (eval(argslist$blupstar)=="cgr")){
    Ztlist <- getME(model,"Ztlist")
    X <- getME(model,"X")
    beta <- getME(model,"beta")

  }

  Ztlist <- c(Ztlist,error=Diagonal(length(blup$ehat)))
  # print(lapply(Ztlist,dim))
  # print(lapply(reff,dim))

  # print(length(Ztlist))
  # print(length(reff))
  # print(length(PBSlist$ranef))

  estar <- mapply(function(pbs,rfi,zti)as.matrix(t(zti)%*%PBS_perm(as.numeric(rfi),pbs)),pbs=PBSlist$ranef,rfi=reff,z = Ztlist,SIMPLIFY = F)
  estar <- matrix(Reduce("+",estar),ncol=np)
  args <- list(estar = estar, assigni = assigni, model = model, X = X, beta = beta, PBSmat = PBSlist$fixeff)


#return(args)

  model0 <- FUN_p(args)

  if(statistic%in% c("quasiF","quasiF_logp") ){
    statp = FUN_stat(model0, assigni)
  }else{
    statp = sapply(model0,function(mod){
    FUN_stat(mod, assigni)
  })}



  out <- list()
  out$model <- model
  out$model0 <- model0
  out$blup <- blup
  out$statistic_distr <- statp
  out$mc <- mc
  out$argslist  <- argslist
  out$PBSlist <- PBSlist
  out$estar <- estar
  out

}




#' Compute quasi F from a from gANOVA_lFormula()
#'@export
lmerModperm.list <- function(model, blupstar = "cgr", np = 4000, method = "terBraak", assigni = 1, statistic = "quasiF_logp",...){


  if(sum(names(model)== c("fr","X","reTrms","REML","formula","wmsgs" ))!=6){
    stop("the model is not the output of gANOVA_lFormula()")
  }

  if((statistic%in%c("Satterthwaite","Satterthwaite_logp","Satterthwaite","Satterthwaite_p"))||(blupstar%in%c("cgr","blup"))){
    cl = quote(gANOVA())
    cl$formula = eval(model$formula)
    cl$data = eval(model$fr)

    model <- eval(cl)
    return(lmerModperm( model = model , blupstar = blupstar, np = np, method = method, assigni = assigni, statistic = statistic,...))


  }


  argslist <- formals()
  mc <- match.call()
  argslist[names(as.list(mc[-1]))] = as.list(mc[-1])


  switch(statistic,
         "quasiF" = {FUN_stat =function(model,assigni){model[,1]}},
         "quasiF_logp" = {FUN_stat =function(model,assigni){abs(pf(q = model[,1] , df1 = model[,10], df2 = model[,11],lower.tail = T, log.p = T))}}
  )

  if(method=="dekker"){stop("the dekker method do not work with quasi F statistics.")}
  FUN_p = function(x)quasiF_terBraak(x)


  blup <- blup_lm(model)

  reff <- blup$Uhat_list
  reff$error <- matrix(blup$ehat,ncol=1)



  PBSlist <- list(
    fixeff = PBSmat(n = nrow(getX(model)), np = np, type = "P"),
    ranef = lapply(lapply(reff,length),function(n){
      PBSmat(n=n,np =np,type = "PS")
    }))


  formula = as.formula(model$formula)
  formula_aov = lme4formula2aov(formula)
  Ztlist <-getZtlm(formula_aov,"reduced",data=model$fr)
  X <- getX(model)
  beta <- qr.coef(qr(cbind(X,t(do.call("rbind",Ztlist)))),model$fr[,attr(terms(model$fr),"response")])[1:ncol(X)]
  Ztlist <- getZtlm(formula_aov,"full",data=model$fr)
  Ztlist <- c(Ztlist,error=Diagonal(length(blup$ehat)))


  estar <- mapply(function(pbs,rfi,zti)as.matrix(t(zti)%*%PBS_perm(as.numeric(rfi),pbs)),pbs=PBSlist$ranef,rfi=reff,z = Ztlist,SIMPLIFY = F)
  estar <- matrix(Reduce("+",estar),ncol=np)
  args <- list(estar = estar, assigni = assigni, model = model, X = X, beta = beta, PBSmat = PBSlist$fixeff)


  model0 <- FUN_p(args)

  statp <- FUN_stat(model0, assigni)

  out <- list()
  out$model <- model
  out$model0 <- model0
  out$blup <- blup
  out$statistic_distr <- statp
  out$mc <- mc
  out$argslist  <- argslist
  out$PBSlist <- PBSlist
  out$estar <- estar
  out


}

