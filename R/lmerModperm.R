#' Permutation test on lmerMod object
#'
#' @description permutation test on lmerMod object
#'
#' @param model a lmerMod object
#' @param blupstar a character string indicating the type of random effect to permute. Default is "cgr" and "blup" and "lm"
#' @param np an integer indicating the number of permutation/bootstrap/suffle.
#' @param method a character string indicating the permutation/bootstrap/suffle method, default is "terBraak".
#' @param assigni a integer indicating the effect to test. Default is 1.
#' @param statistic a character string indicating the test statistic. default is "satterwhaite".
#' @param ... other arguments
#' @importFrom stats pf anova
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



  blup <- blup_FUN(model)


  switch(method,
         "terBraak" = {FUN_p = function(x)lmerModperm_terBraak(x)},
         "dekker" = {FUN_p = function(x)lmerModperm_dekker(x)})

  switch(statistic,
         "Satterthwaite" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,5]}},
         "Satterthwaite_p" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,6]}},
         "Satterthwaite_logp" = {FUN_stat = function(model,assigni){
           ano = anova(model,type=3,ddf="Satterthwaite")
           abs(pf(q=ano[assigni,5], df1 = ano[assigni,3], df2 = ano[assigni,4], lower.tail = F, log.p = T))
         }})




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
  estar <- Reduce("+",estar)
  args <- list(estar = estar, assigni = assigni, model = model, X = X, beta = beta, PBSmat = PBSlist$fixeff)


#return(args)
  model0 = FUN_p(args)
  statp = NULL
  statp = sapply(model0,function(mod){
    FUN_stat(mod, assigni)
  })



  out <- list()
  out$model <- model
  out$model0 <- model0
  out$blup <- blup
  out$statp <- statp
  out$mc <- mc
  out$argslist  <- argslist
  out

}

lmerModperm.gANOVA_lFormula <- function(model, blup_FUN = blup_cgr, np = 4000, method = "terBraak", assigni = 1, statistic = "Satterthwaite",...){

}

