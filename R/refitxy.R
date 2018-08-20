#'Refit both X and Y of a lme4 model
#'
#'@description Refit both X and Y of a gANOVA model.
#'
#'@param object a gANOVA model
#'@param newresp a response vector
#'@param newx a new X matrix of fixed effect.
#'@param rename.response not used.
#'@param ... other arguments.
#'@export
#'@importFrom lme4 checkConv lmerControl mkLmerDevfun optimizeLmer refit
#'@importFrom methods as
refitxy = function(object, newresp = NULL,newx = NULL, rename.response = FALSE, ...){UseMethod("refitxy")}

refitxy.lmerModgANOVA = function(object, newresp = NULL,newx = NULL, rename.response = FALSE, ...){
  if(is.null(newx)){
    return(refit(object = object, newresp = newresp, rename.response = rename.response, ... = ...))
  }
  ## same control than refit refit
  l... <- list(...)
  ctrl.arg <- NULL
  if ("control" %in% names(l...))
    ctrl.arg <- l...$control
  if (!all(names(l...) %in% c("control", "verbose"))) {
    warning("additional arguments to refit.merMod ignored")
  }
  newrespSub <- substitute(newresp)
  if (is.list(newresp)) {
    if (length(newresp) == 1) {
      na.action <- attr(newresp, "na.action")
      newresp <- newresp[[1]]
      attr(newresp, "na.action") <- na.action
    }
    else {
      stop("refit not implemented for 'newresp' lists of length > 1: ",
           "consider ", sQuote("lapply(object,refit)"))
    }
  }
  control <- if (!is.null(ctrl.arg)) {
    if (length(ctrl.arg$optCtrl) == 0) {
      obj.control <- object@optinfo$control
      ignore.pars <- c("xst", "xt")
      if (any(ign <- names(obj.control) %in% ignore.pars))
        obj.control <- obj.control[!ign]
      ctrl.arg$optCtrl <- obj.control
    }
    ctrl.arg
  }else lmerControl()






  if(!is.null(newresp)){
    rcol = attr(attr(model.frame(object), "terms"), "response")
    object@frame[,rcol ] <- newresp
  }

  lfcall <- object@call
  lfcall$formula <- eval(lfcall$formula)
  lfcall[[1]] <- quote(gANOVA::gANOVA_lFormula)
  lfcall$data <- quote(object@frame)

  ##
  start <- as.numeric(getME(object, "theta"))
  lfcall$start <- start

  glmod <- eval(lfcall)
  ## update assign new x
  if(is.null(attr(newx,"assign"))){
    attr(newx,"assign") <- attr(glmod$X,"assign")
  }

  ## change X
  glmod$X <- newx


  verbose = 0L


  devfun <- do.call(mkLmerDevfun, c(glmod, list(start = start,
                                               verbose = verbose, control = control)))

  opt <- optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
                      boundary.tol = control$boundary.tol, control = control$optCtrl,
                      verbose = verbose, start = start, calc.derivs = control$calc.derivs,
                      use.last.params = control$use.last.params)
  cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
                  lbound = environment(devfun)$lower)
  model <- lme4::mkMerMod(environment(devfun), opt, glmod$reTrms,
                         fr = glmod$fr, mc = lfcall, lme4conv = cc)
  res <- lmerTest:::as_lmerModLT(model, devfun)
  res@call <- lfcall
  res <- as(res, c("lmerModgANOVA"))
  res


}
