#' @importFrom gANOVA refit.lmerModgANOVA
#' @importFrom utils setTxtProgressBar txtProgressBar
lmerModperm_dekker = function(args){
  XD = args$X

  D = XD[,attr(XD,"assign")!=args$assigni,drop=F]

  rX = qr.resid(qr(D),XD[,attr(XD,"assign")==args$assigni,drop=F])
  order_xd = order(c(attr(XD,"assign")[attr(XD,"assign")==args$assigni],
                     attr(XD,"assign")[attr(XD,"assign")!=args$assigni]))

  ##ystar
  ystar = as.numeric(XD%*%getME(args$model,"beta")) + args$estar

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
  mi = args$model
  start = getME(mi,"theta")

  ## progress bar
  pb = txtProgressBar(min=1,max = ncol(ystar),initial = 0,style=3)

  for(i in 1:ncol(ystar)){
    setTxtProgressBar(pb,i)



    prxi = PBS_perm(rX,getpbs(args$PBSmat,i))[[1]]
    xdi = cbind(prxi,D)[,order_xd,drop=F]



    model0[[i]] = refitxy(mi,newresp = ystar[,i],newx = xdi,start = start)
    start = (start*(i) + getME(model0[[i]],"theta")*1)/(i+1)
    #model0[[i]] = eval(cl)
  }
  model0[[1]] = args$model


  model0

}
