rm(list=ls())

library(evd)
library(qgam)
library(rmutil)
library(geometry)
library(scales)
library(geometricMVE)
library(mvtnorm)

par(pty="s")
setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/quantiles_emp/")

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

log.strongdep.res = NULL
log.weakdep.res = NULL
gauss.res = NULL
invlog.res = NULL

nrep=200
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  # log - strong dep
  set.seed(k)
  x<-rbvevd(5000,dep=0.4,mar1=c(0,1,0))
  x<-qexp(evd::pgumbel(x))
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  qr = geometricMVE::QR.2d(r,w,tau=tau,method="empirical")
  log.strongdep.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(log.strongdep.res,file="log_strongdep_quantiles.Rdata")}
  
  # log - weak dep
  set.seed(k)
  x<-rbvevd(5000,dep=0.8,mar1=c(0,1,0))
  x<-qexp(evd::pgumbel(x))
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  qr = geometricMVE::QR.2d(r,w,tau=tau,method="empirical")
  log.weakdep.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(log.weakdep.res,file="log_weakdep_quantiles.Rdata")}
  
  # Gaussian
  set.seed(k)
  x<-qexp(pnorm(rmvnorm(5000,sigma=matrix(c(1,0.8,0.8,1),2,2))))
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  qr = geometricMVE::QR.2d(r,w,tau=tau,method="empirical")
  gauss.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(gauss.res,file="gauss_quantiles.Rdata")}
  
  # inverted logistic
  set.seed(k)
  x<-rbvevd(5000,dep=0.7,mar1=c(1,1,1))
  x<-1/x
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  qr = geometricMVE::QR.2d(r,w,tau=tau,method="empirical")
  invlog.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(invlog.res,file="invlog_quantiles.Rdata")}
}
