#model_frame = mf;formula_f = fformula; wformula = formula_within[[2]];idformula = formula_id[[2]]

mm_z <- function(model_frame, formula_f, wformula, idformula){
  tf = terms(formula_f)
  tw = terms(wformula)
  tid = terms(idformula)

  term.labels_between =  (attr(tf,"term.labels")[attr(tf,"order")==1])[
    !(attr(tf,"term.labels")[attr(tf,"order")==1])%in%(attr(tw,"term.labels")[attr(tw,"order")==1])]

  check_between = !sapply(term.labels_between , function(tlbi){
  sum(apply(xtabs(paste(deparse(idformula),tlbi,sep="+"),model_frame),1,function(rowi){
    sum(rowi==0)>1
  }))==0})
  term.labels_between = term.labels_between[check_between]

  formula_between = formula(paste("~ 0",paste(term.labels_between, collapse = ":"),collapse = "+",sep="+"))
  idformula = update.formula(old = idformula,new = ~ 0 + .)

  mm0_id = model.matrix(idformula,model_frame)
  mm0_b = model.matrix(formula_between,model_frame)
  #debug
  # z = lapply(1:NCOL(mm0_b),function(i){
  #   which_row = which(mm0_b[,i] == 1)
  #   which_col = which(colSums(mm0_id[which_row,])!=0)
  #
  #   z0i = mm0_id[which_row,which_col]%*%contr.sum(length(which_col))
  #   zi = matrix(0, ncol= NCOL(z0i),nrow= NROW(mm0_id))
  #   zi[which_row,] = z0i
  #   attr(zi,"contrasts") = contr.sum(length(which_col))
  #   zi
  # })
  #
  # c = as.matrix(bdiag(lapply(z,function(zi)attr(zi,"contrasts"))))
  # z = do.call("cbind",z)
  # attr(z,"contrasts") = c

  mm0_b<<-mm0_b
  sq_i = seq_len(dim(mm0_b)[2])
  if(length(sq_i)==0){sq_i=0}


  z = lapply(sq_i,function(i){
    if(i==0){which_row = 1:nrow(mm0_b)}else{which_row = which(mm0_b[,i] == 1)}

    which_col = which(colSums(mm0_id[which_row,])!=0)

    z0i = mm0_id[which_row,which_col]%*%contr.sum(length(which_col))
    zi = matrix(0, ncol= NCOL(z0i),nrow= NROW(mm0_id))
    zi[which_row,] = z0i
    contr = matrix(0,nrow = ncol(mm0_id),ncol = length(which_col)-1)
    contr[which_col,] = contr.sum(length(which_col))
    attr(zi,"contrasts") = contr
    zi
  })

  lapply(z,function(zi)attr(zi,"contrasts"))

  c = do.call("cbind",lapply(z,function(zi)attr(zi,"contrasts")))
  z = do.call("cbind",z)
  attr(z,"contrasts") = c
  attr(z,"label") = attr(terms(idformula),"term.labels")

  return(z)
}
