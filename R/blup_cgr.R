#' Compute the CGR residuals from a lmerMod object
#'
#' @description Compute the CGR residuals from a lmerMod object.
#'
#' @param model a lmerMod object.
#' @import Matrix
#' @importFrom lme4 VarCorr ranef getME
#' @importFrom stats resid sd
#' @importFrom gANOVA ranef.lmerModgANOVA
#' @export
blup_cgr <- function(model){
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
