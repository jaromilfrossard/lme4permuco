getInteraction = function(var){
  ftable = attr(terms(formula(paste("~",paste(var,collapse="*"),sep=""))),"factors")

  varlist = list()
  for(i in 1:ncol(ftable)){
    sel = ftable[,i]
    if(sum(sel)==1){
      varlist[[i]] = rownames(ftable)[sel==1]}else{
        varlist[[i]] =paste("interaction(",paste(rownames(ftable)[sel==1],collapse=","),")",sep="")
      }

  }
  names(varlist) =colnames(ftable)
  varlist

}
