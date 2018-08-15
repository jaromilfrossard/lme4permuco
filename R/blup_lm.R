#' Compute the blup as fixed effect of a lmerMod object
#'
#' @description Compute the blup as fixed effect of a lmerMod object
#'
#' @param model a lmerMod object.
#' @export
blup_lm = function(model){
  # # Z = t(getME(model,"Zt"))
  # # X = getME(model,"X")
  # faov = lme4formula2aov(model)
  # fformula = lme4:::getFixedFormula(eval(model@call$formula))
  # fformula_design = delete.response(terms(update.formula(fformula, ~.)))
  #
  # terms <- terms(faov, special = "Error", data = model@frame)
  # ind_error <- attr(terms, "specials")$Error
  # error_term <- lapply(ind_error,function(ie)attr(terms, "variables")[[1 + ie]])
  # error_names <- sapply(error_term,function(et){
  #   len_et <- length(et[[2]])
  #   if(len_et == 3){deparse(et[[2]][[2]])}
  #   else if(len_et == 3){deparse(et[[2]])}
  # })
  # names(error_term) <- error_names
  # # extract within part
  # formula_within <- lapply(error_term, function(et){
  #   len_et <- length(et[[2]])
  #   if(len_et==3){formula(paste("~", deparse(et[[2]][[3]]), collapse = ""),env=NULL,showEnv=F)}
  #   else if(len_et==1){
  #     formula(~ 1,env=NULL,showEnv=F) }})
  # # extract id part
  # formula_id <- lapply(error_term,function(et){
  #   len_et <- length(et[[2]])
  #   if(len_et==3){formula(paste("~", deparse(et[[2]][[2]]), collapse = ""),env=NULL,showEnv=F)}
  #   else if(len_et==1){ formula(paste("~", deparse(et[[2]]), collapse = ""),env=NULL,showEnv=F)}})
  # #all fixed formula
  # afformula <- paste(deparse(fformula_design),
  #                    sapply(error_term,function(et)deparse(et[[2]])),collapse="+",sep="+")
  #
  # # create model frame with contrasts
  # mf_design <- permuco:::changeContrast(model@frame, contr = "contr.sum")
  #
  # mf_id <- lapply(formula_id,function(fid){
  #                  model.frame(formula = fid, data = as.data.frame(lapply(mf_design,
  #                               function(col) {col = as.factor(col)
  #                                              contrasts(col) = contr.sum
  #                                              col})))})
  #
  # mf <- model.frame(formula = afformula, data = mf_design)
  # mf_f <- model.frame(formula = fformula, data = mf_design)
  #
  # ## model frame full
  # mf0 = mf
  # for(i in 1:NCOL(mf0)){
  #   if(is.factor(mf0[,i])){
  #     contrasts(mf0[,i],how.many = length(levels(mf0[,i]))) = contrasts(mf0[,i],contrasts=F)
  #   }else{
  #     mf0[,i]<-mf0[,i]
  #   }
  #
  # }
  #
  #
  # # create model.matrix
  #
  # mm <- sparse.model.matrix(fformula_design, data = mf_f)
  # #mm_id <-mapply(function(fid,fw){mm_z(mf,fformula,fw,fid)},fid=formula_id,fw=formula_within,SIMPLIFY = F)
  # contr_id <-mapply(function(fid,fw){contrast_id(mf,fformula,fw,fid)},fid=formula_id,fw=formula_within,SIMPLIFY = F)
  #
  #
  # # create model.matrix ( full)
  # mm0 = sparse.model.matrix(fformula_design,mf0)
  # mm0_id = lapply(formula_id,function(fid) sparse.model.matrix(fid,mf0)[,-1,drop=F])
  #
  #
  #
  # tf = delete.response(terms(update.formula(fformula, ~.)))
  #
  # link = permucoQuasiF:::link(fformula,formula_within[[1]],formula_within[[2]])
  #
  #
  # which_within = list(sort(unique(link[3, ])),
  #                     sort(unique(link[4, ])))
  # which_within = lapply(which_within, function(ww){
  #   names(ww)= c("intercept",colnames(link))[ww+1]
  #   ww
  # })
  # names(which_within) = names(mm0_id)
  #
  # #add third interaction
  # mm0_id[[paste(names(mm0_id),collapse = ":")]] = t(KhatriRao(t(mm0_id[[1]]),t(mm0_id[[2]]),make.dimnames = T))
  # contr_id[[paste(names(contr_id),collapse = ":")]] = kronecker(contr_id[[1]], contr_id[[2]], make.dimnames = T)
  # which_within[[paste(names(which_within),collapse = ":")]]  = (which_within[[1]])[(which_within[[1]])%in%(which_within[[2]])]
  # names(which_within[[3]]) = c("intercept",colnames(link))[which_within[[3]]+1]
  #
  #
  #
  # Ztlist = list()
  # for (i in 1:length(mm0_id)){
  #   Ztlist[[i]] <- Zmat(mm = mm, mm_id= mm0_id[[i]]%*%contr_id[[i]], which_within =which_within[[i]])
  #   names(Ztlist[[i]]) = paste(names(mm0_id)[i],names(Ztlist[[i]]),sep=":")
  #   }
  # Ztlist = unlist(Ztlist,recursive = F)
  #
  # Zmat(mm = mm, mm_id= mm0_id[[1]]%*%contr_id[[1]], which_within =which_within[[1]])
  #
  X = getME(model,"X")
  #Z0tlist = getZtlistlm(model,"full")
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


  # Ztlist = lapply(zm$zs,function(zi){t(Matrix(zi$z0,sparse = T))})
  # assign = unlist(sapply(1:length(Ztlist), function(i)rep(i,nrow(Ztlist[[i]]))))
  # gamma = lapply(1:length(Ztlist), function(i)matrix(gamma[assign==i],ncol=1))
  # names(gamma) = names(Ztlist)
  # if(sum(sapply(zm$zs,function(zi)dim(zi$z))[2,])==nrow(mm_f)){
  #   out = list(Uhat_list = gamma[-length(gamma)], ehat = gamma[length(gamma)],Ztlist = Ztlist, X = mm_f, beta = beta[-c(1:ncol(mm_f))])
  #   }



}

