#' Permutation test on lmerMod object
#'
#' @description permutation test on lmerMod object
#'
#' @param model a lmerMod object
#' @param blup_FUN the type of BLUP to permute. Default is CGR.
#' @param np an integer indicating the number of permutation/bootstrap/suffle.
#' @param method a character string indicating the permutation/bootstrap/suffle method, default is "terBraak".
#' @param assigni a integer indicating the effect to test. Default is 1.
#' @param statistic a character string indicating the test statistic. default is "satterwhaite".
#' @param ... other arguments
#' @importFrom stats pf anova
#' @export
lmerModperm <- function(model, blup_FUN, np, method, assigni, statistic, ...){UseMethod("lmerModperm")}

#' @export
lmerModperm.lmerModgANOVA <- function(model, blup_FUN = blup_cgr, np = 4000, method = "terBraak", assigni = 1, statistic = "Satterthwaite",...){



  blup <- blup_FUN(model)


  switch(method,
         "terBraak" = {FUN_p = function(x)lmerModperm_terBraak(x)},
         "dekker" = {FUN_p = function(x)lmerModperm_terBraak(x)})

  switch(statistic,
         "Satterthwaite" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,5]}},
         "Satterthwaite_p" = {FUN_stat = function(model,assigni){anova(model,type=3,ddf="Satterthwaite")[assigni,6]}},
         "Satterthwaite_logp" = {FUN_stat = function(model,assigni){
           ano = anova(model,type=3,ddf="Satterthwaite")
           abs(pf(q=ano[assigni,5], df1 = ano[assigni,3], df2 = ano[assigni,4], lower.tail = F, log.p = T))
         }})




  ### create new error
  Ztlist <- getME(model,"Ztlist")

  reff <- blup$Uhat_list
  reff$error <- matrix(blup$ehat,ncol=1)
  Ztlist <- c(Ztlist,error=Diagonal(length(blup$ehat)))

  PBSlist <- lapply(lapply(reff,length),function(n){
    PBSmat(n=n,np =np,type = "PS")
  })



  estar <- mapply(function(pbs,rfi,zti)as.matrix(t(zti)%*%PBS_perm(as.numeric(rfi),pbs)),pbs=PBSlist,rfi=reff,z = Ztlist,SIMPLIFY = F)
  estar <- Reduce("+",estar)
  fitted <- fitted (model)
  args <- list(estar = estar, fitted = fitted, assigni = assigni, model = model)


#return(args)
  model0 = FUN_p(args)
  statp = NULL
  statp = sapply(model0,function(mod){
    FUN_stat(mod, assigni)
  })



  out <- list()
  out$model <- model
  out$model0 = model0
  out$blup <- blup
  out$statp = statp
  out

}
