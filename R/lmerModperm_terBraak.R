#' @importFrom gANOVA refit.lmerModgANOVA
#'
lmerModperm_terBraak = function(args){
  XD = args$X
  beta = args$beta

  fitted_star = as.numeric(XD[,attr(XD,"assign")!=args$assigni,drop=F]%*%beta[attr(XD,"assign")!=args$assigni])


  ## fitted value
  ystar=args$estar+fitted_star
  ##sarting valeu
  theta = getME(args$model,"theta")
  cl = args$model@call
  cl[[1]]=gANOVA::gANOVA
  cl$start = eval(theta)

  ## new parameter
  cl$formula = eval(cl$formula)
  f = eval(cl$formula)
  model0 = list()
  length(model0)= ncol(ystar)

  ## progress bar
  pb = txtProgressBar(min=1,max = ncol(ystar),initial = 0,style=3)


  for(i in 1:ncol(ystar)){
    setTxtProgressBar(pb,i)
    model0[[i]] = refit.lmerModgANOVA(args$model,newresp = ystar[,i], start = start)
    start = (start*(i) + getME(model0[[i]],"theta")*1)/(i+1)

  }
  model0[[1]] = args$model


  model0

}
