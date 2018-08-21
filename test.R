rm(list=ls())
library(lme4)
library(lmerTest)
library(MASS)
library(combinat)
library(gANOVA)
#devtools::install_github("jaromilfrossard/gANOVA")
Sys.setenv(LANG = "en")

dir= ""

#### simulation type
## decreasing variance through interaction
## increasing variable through interaction
## Large variance for id:item
##



t0 = proc.time()
#setwd(dir)


lf= list.files(paste(dir,"function_data",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"function_data/",lfi,sep=""))
}


ni= 18 #18
ns = 18#16
na =list(S = 3,I =3, C =3, SI = 1, O = 1, R = 1)


fseed = 42

source(paste("data_uni_lmm.R",sep=""))




f_s = subject~(fI*fC)
f_i = item~(fS*fC)
f_si = (subject_item)~(fC)

zl_s = Zlist(f_s,df)
scale_s = scalelist(Zlist = zl_s,type="flat")
omega_glist_s = Omega_glist(Zlist = zl_s,scale = scale_s, type = "homoscedastic")


zl_i = Zlist(f_i,df)
scale_i = scalelist(Zlist = zl_i,type="flat")
omega_glist_i = Omega_glist(Zlist = zl_i,scale = scale_i, type = "homoscedastic")


if(!is.null(f_si)){
  zl_si = Zlist(f_si,df)
  scale_si = scalelist(Zlist = zl_si,type="decreasing")
  omega_glist_si = Omega_glist(Zlist = zl_si,scale = scale_i, type = "homoscedastic")
  Omega = omega_glist_s$Omega + omega_glist_i$Omega + omega_glist_si$Omega+
    Diagonal(n = nrow(df))*(max(scale_s,scale_i,scale_si)*1.5)^2}else{
      Omega = omega_glist_s$Omega + omega_glist_i$Omega +
        Diagonal(n = nrow(df))*(max(scale_s,scale_i)*1.5)^2}

err = MASS:::mvrnorm(n=1, mu= rep(0,NCOL(Omega)), Sigma=Omega)
df$y=err

#matrixcalc:::is.positive.definite(as.matrix(Omega))


df =permuco:::changeContrast(df,contr = contr.sum)

df=cbind(df,model.matrix(~fC+fI+fS,df,contrasts.arg = list(fC = contr.poly, fS = contr.poly, fI=contr.poly))[,-1])




ffull = y~fS*fC*fI + (1|subject)+ (1|subject:fC)+ (1|subject:fI) +
  (1|subject:fC:fI) +
  (1|item) + (1|item:fC)+ (1|item:fS)+ (1|item:fC:fS)+
  (1| subject_item)

fsmall = y~fS*fC*fI + (1|subject)+ (1|subject:fC)+ (1|item)
fmedium= y~fS*fC*fI  + (1|subject)+ (1|subject:fC)+ (1|subject:fI) +
  (1|subject:fC:fI)+ (1|item)+ (1|item:fC)
flarge= y~fS*fC*fI  + (1|subject)+ (1|subject:fC)+ (1|subject:fI) +
  (1|subject:fC:fI)+ (1|item)+ (1|item:fC)+ (1|item:fS)+ (1|item:fC:fS)
fmlmer= y~fS*fC*fI  + ((fS.L+ fS.Q)||subject)

#devtools::install_github("jaromilfrossard/lme4permuco")
library(lme4permuco)
modelg = gANOVA(fmedium,df,REML =T)
#modelg = gANOVA(flarge,df,REML =T)
# object = modelg
# newresp= rev(getME(modelg,"y"))
# newx = getME(modelg,"X")[sample(length(newresp)),]


# lf= list.files(paste(dir,"R",sep=""))
#
# for(lfi in lf){
#   print(lfi)
#   source(paste(dir,"R/",lfi,sep=""))
# }


params = expand.grid(method = c("dekker","terBraak"),
blup = c("lm","blup","cgr"),stringsAsFactors = F)

modlist = list()

i= 2
np = 200

for(i in 1: nrow(params)){
  modlist[[i]] = lmerModperm(modelg,np=np,blupstar = params$blup[i],method = params$method[i])}


perm_tb_lm = lmerModperm(modelg,np=10,blup_FUN = blup_lm)
perm_tb_blup = lmerModperm(modelg,np=200,blup_FUN = blup_blup)
perm_tb_cgr = lmerModperm(modelg,np=200,blup_FUN = blup_cgr)

perm_tb_lm = lmerModperm(modelg,np=200,blup_FUN = blup_lm)
perm_tb_blup = lmerModperm(modelg,np=200,blup_FUN = blup_blup)
perm_tb_cgr = lmerModperm(modelg,np=200,blup_FUN = blup_cgr)


getSD = function(model){
  vc = VarCorr(model)
  sapply(vc, function(vci){
    attr(vci, "stddev")
  })
}

##### plot thetas
thetas_lm = t(sapply(permlm$model0,function(mod){getSD(mod)}))
thetas_cgr = t(sapply(permcgr$model0,function(mod){getSD(mod)}))
thetas_blup = t(sapply(permblup$model0,function(mod){getSD(mod)}))


sapply(permlm$model0,function(mod)mod@optinfo$conv$opt==0)
sapply(permcgr$model0,function(mod)mod@optinfo$conv$opt==0)
sapply(permblup$model0,function(mod)mod@optinfo$conv$opt==0)


ntheta = ncol(thetas_lm)
par(mfcol=c(3,ntheta),oma=c(0,0,0,0),mar=c(2,2,1,.5))

for(i in 1:ntheta){
  xlim= range(c(thetas_lm[,i],thetas_cgr[,i],thetas_blup[,i]))
plot(density(thetas_lm[,i]),main="",xlab = "",ylab="",xlim=xlim)
abline(v= thetas_lm[1,i])
plot(density(thetas_cgr[,i]),main="",xlab = "",ylab="",xlim=xlim)
abline(v= thetas_cgr[1,i])
plot(density(thetas_blup[,i]),main="",xlab = "",ylab="",xlim=xlim)
abline(v= thetas_blup[1,i])}



a$model

i=2
theta = getME(a$model,"theta")
cl = a$model@call
cl[[1]]=gANOVA::gANOVA
cl$start = eval(theta)
cl$formula = eval(cl$formula)
f = eval(cl$formula)
f = formula(paste("ystar[,",eval(i),"]",as.character(f)[1],as.character(f)[3]))
cl$formula = f
model0 = eval(cl)
gANOVA:::anova.lmerModgANOVA(model0,"satterwhaite")

lmerModperm_terBraak(args)



