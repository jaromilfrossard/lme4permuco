rm(list=ls())
library(lme4)
library(lmerTest)
library(MASS)
library(combinat)
library(gANOVA)
devtools::install_github("jaromilfrossard/permucoQuasiF")
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
na =list(S = 3,I =3, C =2, SI = 1, O = 1, R = 1)


fseed = 42

source(paste("data_uni_lmm.R",sep=""))


head(df)

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
fqf=y ~ fS * fC * fI + Error(subject/(fI*fC)) + Error(item/(fS*fC))

fsmall = y~fS*fC*fI + (1|subject)+ (1|subject:fC)+ (1|item)
fmedium= y~fS*fC*fI  + (1|subject)+ (1|subject:fC)+ (1|subject:fI) +
  (1|subject:fC:fI)+ (1|item)+ (1|item:fC)
flarge= y~fS*fC*fI  + (1|subject)+ (1|subject:fC)+ (1|subject:fI) +
  (1|subject:fC:fI)+ (1|item)+ (1|item:fC)+ (1|item:fS)+ (1|item:fC:fS)
fmlmer= y~fS*fC*fI  + ((fS.L+ fS.Q)||subject)


f_g_si = formula(paste("y~fS*fC*fI",
                       "+(1|item) + (1|item:fS)+ (1|item:fC)+ (1|item:fS:fC)",
                       "+(1|subject) + (1|subject:fI)+ (1|subject:fC) + (1|subject:fI:fC)",
                       "+(1|subject_item)",sep=""))




#devtools::install_github("jaromilfrossard/lme4permuco")
#library(lme4permuco)
df2 = droplevels(df[df$fC!="c_c",])
df3 = df[-c(1,50,100,150,200,220),]
m1=gANOVA(fmedium,df3,REML=T)
m2=lmer(fmedium,df3,REML=T)


lf= list.files(paste(dir,"R",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"R/",lfi,sep=""))
}




model = gANOVA_lFormula(f_g_si,df,REML =T)
modelg = gANOVA(f_g_si,df,REML =T)




qf_tb_lm = lmerModperm(modelg,np = 4000,method = "terBraak",statistic = "quasiF_logp",blupstar = "lm")
qf_tb_cgr = lmerModperm(modelg,np = 4000,method = "terBraak",statistic = "quasiF_logp",blupstar = "cgr")
qf_tb_blup = lmerModperm(modelg,np = 4000,method = "terBraak",statistic = "quasiF_logp",blupstar = "blup")

f = lme4formula2aov(f_g_si)
f[[3]] = f[[3]][[2]]

head(qf_tb_lm$model0)

permucoQuasiF::aovperm(f,df)


qf_tb_lm$statistic_distr[1]
qf_tb_cgr$statistic_distr[1]
qf_tb_blup$statistic_distr[1]

head(qf_tb_lm$model0)
head(qf_tb_cgr$model0)

terms(fformula)
formula_within



fformula = formula(paste(deparse(f0[[2]]),"~",deparse(f0[[3]][[2]][[2]][[2]]),sep=""))


sapply(Ztlist,dim)

qf$model

(proc.time()-t0)/100

qf = lmerModperm(modelg,np=20)

thetas = t(sapply(qf$model0,function(mod){getME(mod,"theta")}))
cbind(apply(t(t(thetas)-thetas[1,]),1,function(x)sum(x^2)),
apply(t(t(thetas)-colMeans(thetas)),1,function(x)sum(x^2)))



qf = lmerModperm(modelg,np=20,method = "dekker")

colMeans(t(sapply(qf$model0,function(mi){getME(mi,"theta")})))

apply(t(sapply(qf$model0,function(mi){getME(mi,"theta")})),2,cumsum)/1:20


lfg = gANOVA_lFormula(fmedium,df,REML =T)

A = lmerModperm(lfg,np=20)




qflp = lmerModperm(modelg,np=20,blupstar = "lm",method = "terBraak",statistic = "quasiF_logp")


paramSat = expand.grid(method = c("dekker","terBraak"),
blup = c("lm","blup","cgr"), statistic= "Satterthwaite",stringsAsFactors = F)
paramqf = expand.grid(method = c("terBraak"),
                      blup = c("lm","blup","cgr"), statistic= c("quasiF","quasiF_logp"),stringsAsFactors = F)

params = rbind(paramSat,paramqf)
modlist = list()


lmerModperm(modelg,np=3,blupstar = params$blup[1],method = params$method[1],statistic =params$statistic[1])

i= 2
np = 200

for(i in 1: nrow(params)){
  print(params[i,])
  modlist[[i]] = lmerModperm(modelg,np=np,blupstar = params$blup[i],method = params$method[i],statistic =params$statistic[i] )}


save(modlist,file="modlist.RData")

getSD = function(model){
  vc = VarCorr(model)
  sapply(vc, function(vci){
    attr(vci, "stddev")
  })
}

##### plot thetas
sds = lapply(modlist[1:6],function(modi)t(sapply(modi$model0,function(mod){getSD(mod)})))


nsd = ncol(sds[[1]])
par(mfcol=c(length(sds),ntheta),oma=c(2,2,0,0),mar=c(2,2,1,.5))

for(sdi in 1:nsd){
  xlim = range(sapply(sds,function(sdmat){range(sdmat[,sdi])}))
  for(modi in 1:length(sds)){
    xi = sds[[modi]][,sdi]
    plot(density(xi),main="",xlab = "",ylab="",xlim=xlim)
    abline(v= xi[1])
  }
}



  xlim= range(c(thetas_lm[,i],thetas_cgr[,i],thetas_blup[,i]))

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



