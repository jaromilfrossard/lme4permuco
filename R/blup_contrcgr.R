#' Compute the contrasts of the CGR residuals from a lmerMod object modified
#'
#' @description Compute the contrasts of the CGR residuals from a lmerMod object.
#'
#' @param model a gANOVA object.
#' @export
blup_contrcgr <- function(model){
  model.ranef <- ranef(model)
  model.resid <- resid(model)
  Uhat.list <- lapply(seq_along(model.ranef), FUN = function(i) {
    u <- scale(model.ranef[[i]], scale = FALSE)
    if(sd(u)<=1e-8){
      return(u)
    }
    S <- (t(u) %*% u)/length(u)
    re.name <- names(model.ranef)[i]
    R <- bdiag(VarCorr(model)[[names(model.ranef)[i]]])
    Ls <- chol(S, pivot = TRUE)
    Lr <- chol(R, pivot = TRUE)
    A <- t(Lr %*% solve(Ls))
    Uhat <- as.matrix(u %*% A)
    return(Uhat)
  })
  names(Uhat.list) <- names(model.ranef)

  ## modify by contratsts
  contrlist <- getContrlist(model)

  #length of sampling unit for full contrasts
  SUn <- lapply(names(model@reTrms$contrlist),function(namei) {
    namei = unlist(strsplit(unlist(strsplit(namei, "[|]"))[2],"[:]"))[1]
    length(levels(droplevels(model@frame[,gsub(" ", "", namei, fixed = TRUE)])))

  })
  Uhat.list <- mapply(function(ni,contri,uhati)(Diagonal(ni)%x%contri)%*%uhati,
                      contri=contrlist,ni=SUn, uhati = Uhat.list,SIMPLIFY = F)



  e <- as.numeric(scale(model.resid, scale = FALSE))
  sigma <- lme4::getME(model, "sigma")
  ehat <- sigma * e * (sum(e^2)/length(e))^(-1/2)
  #Z <- lme4::getME(object = model, name = "Ztlist")
  #X <- getME(model,"X")
  #beta <- getME(model,"beta")
  #level.num <- lme4::getME(object = model, name = "n_rfacs")

  #out = list(X = X, beta = beta, Z = Z, Uhat_list = Uhat.list, ehat = ehat)
  out = list(Uhat_list = Uhat.list, ehat = ehat)
  return(out)

}
