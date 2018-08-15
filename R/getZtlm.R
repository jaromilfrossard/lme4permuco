#' Extract / create the Z matrix of type LM from a model
#'
#'@description Extract / create the Z matrix of type LM from a model
#'@param model a lmerModgANOVA model
#'@param type integer indicating the type of matrix "full" for Z0, "reduced" for Z = Z0C' or "contrasts" for C
#'@export
#'@importFrom stats delete.response model.frame contrasts 'contrasts<-'
getZtlm <- function(model, type = "full"){UseMethod("getZtlm")}

#'@export
getZtlm <- function(model, type){
  # Z = t(getME(model,"Zt"))
  # X = getME(model,"X")
  faov = lme4formula2aov(model)
  fformula = lme4:::getFixedFormula(eval(model@call$formula))
  fformula_design = delete.response(terms(update.formula(fformula, ~.)))

  terms <- terms(faov, special = "Error", data = model@frame)
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
  # extract id part
  formula_id <- lapply(error_term,function(et){
    len_et <- length(et[[2]])
    if(len_et==3){formula(paste("~", deparse(et[[2]][[2]]), collapse = ""),env=NULL,showEnv=F)}
    else if(len_et==1){ formula(paste("~", deparse(et[[2]]), collapse = ""),env=NULL,showEnv=F)}})
  #all fixed formula
  afformula <- paste(deparse(fformula_design),
                     sapply(error_term,function(et)deparse(et[[2]])),collapse="+",sep="+")

  # create model frame with contrasts
  mf_design <- permuco:::changeContrast(model@frame, contr = "contr.sum")

  mf_id <- lapply(formula_id,function(fid){
    model.frame(formula = fid, data = as.data.frame(lapply(mf_design,
                                                           function(col) {col = as.factor(col)
                                                           contrasts(col) = contr.sum
                                                           col})))})

  mf <- model.frame(formula = afformula, data = mf_design)
  mf_f <- model.frame(formula = fformula, data = mf_design)

  ## model frame full
  mf0 = mf
  for(i in 1:NCOL(mf0)){
    if(is.factor(mf0[,i])){
      contrasts(mf0[,i],how.many = length(levels(mf0[,i]))) = contrasts(mf0[,i],contrasts=F)
    }else{
      mf0[,i]<-mf0[,i]
    }

  }



  # create model.matrix

  mm <- sparse.model.matrix(fformula_design, data = mf_f)
  #mm_id <-mapply(function(fid,fw){mm_z(mf,fformula,fw,fid)},fid=formula_id,fw=formula_within,SIMPLIFY = F)
  contr_id <-mapply(function(fid,fw){contrast_id(mf,fformula,fw,fid)},fid=formula_id,fw=formula_within,SIMPLIFY = F)


  # create model.matrix ( full)
  mm0 = sparse.model.matrix(fformula_design,mf0)
  mm0_id = lapply(formula_id,function(fid) sparse.model.matrix(fid,mf0)[,-1,drop=F])



  tf = delete.response(terms(update.formula(fformula, ~.)))

  link = permucoQuasiF:::link(fformula,formula_within[[1]],formula_within[[2]])


  which_within = list(sort(unique(link[3, ])),
                      sort(unique(link[4, ])))
  which_within = lapply(which_within, function(ww){
    names(ww)= c("intercept",colnames(link))[ww+1]
    ww
  })
  names(which_within) = names(mm0_id)

  #add third interaction
  mm0_id[[paste(names(mm0_id),collapse = ":")]] = t(KhatriRao(t(mm0_id[[1]]),t(mm0_id[[2]]),make.dimnames = T))
  contr_id[[paste(names(contr_id),collapse = ":")]] = kronecker(contr_id[[1]], contr_id[[2]], make.dimnames = T)
  which_within[[paste(names(which_within),collapse = ":")]]  = (which_within[[1]])[(which_within[[1]])%in%(which_within[[2]])]
  names(which_within[[3]]) = c("intercept",colnames(link))[which_within[[3]]+1]



  ####fixed contrast

  if(type=="contrasts"){
  contr_mm = list()
  contr_fact = mapply(function(ctype,cdim){
    dn = dimnames(cdim)
    dn[[2]] = dn[[2]][1:(length(dn[[2]])-1)]
    ci = eval(parse(text = ctype))(ncol(cdim))
    dimnames(ci) = dn
    ci},
    ctype = attr(mm,"contrast"),
    cdim = attr(mm0,"contrast"),SIMPLIFY = F)

  for(i in 1:ncol(attr(terms(fformula_design),"factors"))){
    contri = contr_fact[attr(terms(fformula_design),"factors")[,i]]
    contri = lapply(contri,function(ci)Matrix(ci,dimnames = dimnames(ci),sparse = T))
    contr_mm[[i]] = Reducekronecker(contri,make.dimnames = T)
  }
  names(contr_mm) = colnames(attr(terms(fformula_design),"factors"))}




  Ztlist = list()
  switch(type,
         "full" ={
           for (i in 1:length(mm0_id)){
             Ztlist[[i]] <- Ztmat(mm = mm0, mm_id= mm0_id[[i]], which_within =which_within[[i]])
             names(Ztlist[[i]]) = paste(names(mm0_id)[i],names(Ztlist[[i]]),sep=":")
           }
           Ztlist = unlist(Ztlist,recursive = F)},
         "reduced" ={
           for (i in 1:length(mm0_id)){
             Ztlist[[i]] <- Ztmat(mm = mm, mm_id= mm0_id[[i]]%*%contr_id[[i]], which_within =which_within[[i]])
             names(Ztlist[[i]]) = paste(names(mm0_id)[i],names(Ztlist[[i]]),sep=":")
           }
           Ztlist = unlist(Ztlist,recursive = F)},
         "contrasts" = {
           for (i in 1:length(mm0_id)){
             Ztlist[[i]] = lapply(which_within[[i]],function(ww){
               if(ww==0){
                 contr_id[[i]]
               }else{
                 kronecker(contr_id[[i]],contr_mm[[ww]],make.dimnames =T)

               }
             })
             names(Ztlist[[i]]) = paste(names(mm0_id)[i],names(Ztlist[[i]]),sep=":")
           }
           Ztlist = unlist(Ztlist,recursive = F)})
  return(Ztlist)
}
