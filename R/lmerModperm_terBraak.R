#' @importFrom gANOVA refit.lmerModgANOVA
#'
lmerModperm_terBraak = function(args){
  X = getME(args$model,"X")
  beta = getME(args$model,"beta")

  fitted_star = as.numeric(X[,attr(X,"assign")!=args$assigni,drop=F]%*%beta[attr(X,"assign")!=args$assigni])


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

  prog = 0L
  cat( paste (prog, "%\n"))

  for(i in 1:ncol(ystar)){
    if(prog < as.integer(round(i/ncol(ystar)*100))){
      prog = as.integer(round(i/ncol(ystar)*100))
      cat( paste (prog, "%\n"))
    }

  # f = paste("ystar[,",eval(i),"]",as.character(f)[1],as.character(f)[3])
  # cl$formula = update.formula(old= formula(cl$formula),new = eval(parse(text=f)))
  # print(cl)
  model0[[i]] = refit.lmerModgANOVA(args$model,newresp = ystar[,i])
  #model0[[i]] = eval(cl)
  }
  model0[[1]] = args$model


  model0

}
