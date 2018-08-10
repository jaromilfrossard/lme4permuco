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
  (1|subject:fC:fI)
fmlmer= y~fS*fC*fI  + ((fS.L+ fS.Q)||subject)

modelg = gANOVA(fmedium,df,REML =T)
modelmg = lmer(fmlmer,df,REML =T)



model = modelg
frfull=ranefcontrasts(model)



formula <- formula(model)
bars <- findbars(lme4:::RHSForm(formula))
names(bars) <- lme4:::barnames(bars)
fr <-  model@frame
rf <- ranef(model)

bars <- bars[match(names(bars),names(rf))]

contrlist <- lapply(bars, mkcontrastslist, fr)


lapply(contrlist,function(x)dim(x$sm))

i=1
as.numeric(contrlist[[i]]$contr%*%matrix(as.numeric(as.matrix(rf[[i]])),ncol=blist[[i]]$n))
dim(contrlist[[1]]$sm)



sapply(contrlist[[1]]$sm,dim)

rf$`subject:fC:fI`[1:4,]
lapply(blist,function(x)x$n)

a$contr
str(sm)
rf$`subject:fC:fI`

do.call("%x%",a$contr)

a$contr
class(a$sm)


names(a$sm)

image(a$sm[[1]])

lme4:::RHSForm(formula)
i=2
sumrf=t(ztlistl[[1]])%*%(as.matrix(rfl[[1]]))+
t(ztlistl[[2]])%*%(as.matrix(rfl[[2]]))+
t(ztlistl[[3]])%*%(as.matrix(rfl[[3]]))
sumrf-getME(modell,"Z")%*%getME(modell,"b")
cbind(getME(modell,"u"),getME(modell,"b"))
cbind(getME(modelg,"u"),getME(modelg,"b"))

ranef(model)

formula(model)


library(lme4permuco)

x = abs(matrix(rnorm(10*5),nrow=10))
PBSmat = PBSmat(10,12,"PS")
PBS_perm(x = x, PBSmat)

lf= list.files("R")
for(lfi in lf){
  print(lfi)
  source(paste("R/",lfi,sep=""))
}

lmerp = lmerModperm(model)


str(lmerp$blup)

lme4::VarCorr(model)

lme4:::VarCorr.merMod

np = 500
method = "dekker"
blup_FUN =blup_cgr
statistic ="quasif"
assigni = 1
t0 = proc.time()
qf_cgr_d_reml = randomize(model=model,np= np,method = "dekker",blup_FUN =blup_cgr,statistic ="quasif")


