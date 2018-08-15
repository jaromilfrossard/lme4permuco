#mm = mm; mm_id = mm0_id[[1]]%*%contr_id[[1]]; linki =link[3,]

Ztmat <- function(mm, mm_id, which_within ){

  assign = attr(mm, "assign")


  zm = lapply(which_within,function(ww){
    KhatriRao(t(mm_id),t(mm[,assign==ww,drop=F]),make.dimnames = T)
  })
  names(zm) <- names(which_within)
  zm
}
