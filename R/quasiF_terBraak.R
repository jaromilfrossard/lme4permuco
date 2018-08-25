quasiF_terBraak <- function(args){
  if(class(args$model)=="lmerModgANOVA"){
    data = args$model@frame
    formula <- formula(args$model)
  }else if(class(args$model)=="list"){
    data = args$model$fr
    formula <- args$model$formula
  }
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
  formula_full = formula(paste(deparse(formula[[2]]),formula_full))


  #create link
  terms <- terms(formula_full, special = "Error", data = data)
  terms <- delete.response(terms)
  ind_error <- attr(terms, "specials")$Error

  error_term <- lapply(ind_error,function(ie)attr(terms, "variables")[[1 + ie]])
  error_names <- sapply(error_term,function(et){
    len_et <- length(et[[2]])
    if(len_et == 3){deparse(et[[2]][[2]])}
    else if(len_et == 3){deparse(et[[2]])}
  })
  names(error_term) <- error_names
  # extract within part
  formula_within <- lapply(error_term, function(et){
    len_et <- length(et[[2]])
    if(len_et==3){formula(paste("~", deparse(et[[2]][[3]]), collapse = ""),env=NULL,showEnv=F)}
    else if(len_et==1){
      formula(~ 1,env=NULL,showEnv=F) }})


  link = permucoQuasiF:::link(formula_fix,formula_within[[1]],formula_within[[2]])
  #link ok

  ## creating matrix
  formula_su = formula(paste(c(deparse(formula_fix),paste("Error(",SU,"/1)",sep="")),collapse = "+",sep="+"))

  Ztlist_su =getZtlm(formula_su,"reduced",data=data)

  assign = attr(args$X, "assign")
  select_x = assign == args$assigni
  select_w = list(w1 = assign == (link[3, args$assigni]),
                  w2 = assign == (link[4, args$assigni]))
  select_w$select_w3 =select_w$select_w1&select_w$select_w2

  zs = mapply(function(zi,wi){
    if( sum(wi)==0){zi}else{
    KhatriRao(zi,t(args$X[,wi,drop=F]))}
  },zi = Ztlist_su, wi = select_w,SIMPLIFY = F)


  #terbraak Y

  XD = args$X
  beta = args$beta

  fitted_star = as.numeric(XD[,attr(XD,"assign")!=args$assigni,drop=F]%*%beta[attr(XD,"assign")!=args$assigni])


  ## fitted value
  ystar=args$estar+fitted_star
  ystar[,1] = ystar[,1] + as.numeric(XD[,attr(XD,"assign")==args$assigni,drop=F]%*%beta[attr(XD,"assign")==args$assigni])
  ####computing QR


  #qr_d = qr(args$X[, !select_x, drop = F])
  rdx <- Matrix(XD[, select_x, drop = F],sparse=T)
  qr_rdx <- Matrix::qr(rdx)


  qr_z <- lapply(zs,function(zi){
    qr(t(zi))
  })


  qr_all <- c(num1 = qr_rdx,
             num2 = qr_z[[3]],
             den1 = qr_z[[1]],
             den2 = qr_z[[2]])

  rank_all = list(num1 = sum(select_x),
               num2 = nrow(zs[[3]]),
               den1 = nrow(zs[[1]]),
               den2 = nrow(zs[[2]]))




  SS = lapply(qr_all,function(qri){
    Matrix::colSums(qr.fitted(qri, ystar)^2)
  })

  MS = mapply(function(ssi,ri){
    ssi/ri
  },ssi = SS,ri = rank_all,SIMPLIFY = F)



  dfn = (MS$num1+MS$num2)^2/(MS$num1^2/rank_all$num1+MS$num2^2/rank_all$num2)
  dfd = (MS$den1+MS$den2)^2/(MS$den1^2/rank_all$den1+MS$den2^2/rank_all$den2)

  quasif = (MS$num1+MS$num2)/(MS$den1+MS$den2)

  out = cbind(quasif = quasif, SSn1 = SS$num1, SSn2 = SS$num2, SSd1 = SS$den1, SSd2 = SS$den2,
          dfn1 = rank_all$num1,dfn2 = rank_all$num2,dfd1 = rank_all$den1, dfd2 = rank_all$den2,dfn = dfn ,dfd = dfd,
           `parametric P(>F)` = pf(q = quasif , df1 = dfn, df2 = dfd,lower.tail = T))
  out




}

