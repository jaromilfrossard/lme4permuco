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


library(lme4permuco)
modelg = gANOVA(fmedium,df,REML =T)
#modelg = gANOVA(flarge,df,REML =T)


permlm = lmerModperm(modelg,np=200,blup_FUN = blup_lm)
permblup = lmerModperm(modelg,np=200,blup_FUN = blup_blup)
permcgr = lmerModperm(modelg,np=200,blup_FUN = blup_cgr)

lmc= as.list(mc[-1])

bllm=blup_lm(model)

bllm$ehat

lapply(bllm$Uhat_list,dim)

lme4:::getFixedFormula()

library(lme4permuco)

permucoQuasiF:::clusterlm_quasif


lf= list.files(paste(dir,"R",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"R/",lfi,sep=""))
}


lf= list.files(paste(dir,"../gANOVA/R",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"../gANOVA/R/",lfi,sep=""))
}




model = modelg



blup = lmerModperm(modelg,np=200,blup_FUN = blup_blup,statistic = "Satterthwaite")
cgr = lmerModperm(modelg,np=200,blup_FUN = blup_cgr,statistic = "Satterthwaite")



mean(blup$statp>=(blup$statp[1]))
mean(cgr$statp>=(cgr$statp[1]))

cgr$statp[1]
blup$statp[1]
cgr
blupp$statp
bluplp$statp

getME(blup$model0[[2]])
lapply(blup$model0, function(mod)lmerTest:::anova.lmerModLmerTest(mod,ddf = "Satterthwaite",type=3))

a = lmerTest:::anova.lmerModLmerTest(blup$model0[[2]])
lmerTest:::contestMD.lmerModLmerTest

FUN_stat(blup$model0[[2]],2)


getME(blup$model0[[2]],"L")
getME(blup$model0[[2]],"beta")

lmerTest:::as_lmerModLmerTest(blup$model0[[2]])
refit(modelg,newresp = rev(df$y))

getME(model,"theta")
anova(model,type = "III",method = "satterthwaite")

lmerTest:::anova.lmerModLmerTest

blup = lmerModperm(modelg,np=3,blup_FUN = blup_blup)
blup$statp

cgr = lmerModperm(model,np=200)



thetas_cgr = t(sapply(cgr$model0,function(mod){getME(mod,"theta")}))
thetas_blup = t(sapply(blup$model0,function(mod){getME(mod,"theta")}))


model.ranef <- ranef(model)
model.resid <- resid(model)

par(mfrow=c(1,2))
i = 5
plot(density(thetas_cgr[,i]))
abline(v= thetas_cgr[1,i])
plot(density(thetas_blup[,i]))
abline(v= thetas_blup[1,i])

sapply(a$statp,function(mod){getME(mod,"theta")})



a$

ystar=a$estar+a$fitted


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



