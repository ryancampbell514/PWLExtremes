rm(list=ls())

library(evd)
# library(qgam)
library(rmutil)
library(geometry)
library(scales)
# library(geometricMVE)
library(mvtnorm)
library(this.path)

# par(pty="s")
# setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/quantiles/")
# 
# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# library(this.path)
# setwd(this.path::here())
setwd(this.path::here())

fn.dir = "../../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# source("../../../kde_quants_bw_plots_unconstrained/thresh_estimation_unconstrained.R")
source("../../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")

log.strongdep.res = list(list(NULL))
log.weakdep.res = list(list(NULL))
gauss.res = list(list(NULL))
invlog.res = list(list(NULL))

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
  # qr = QR.2d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  log.strongdep.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(log.strongdep.res,file="log_strongdep_quantiles.Rdata")}
  
  # log - weak dep
  set.seed(k)
  x<-rbvevd(5000,dep=0.8,mar1=c(0,1,0))
  x<-qexp(evd::pgumbel(x))
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  # qr = QR.2d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  log.weakdep.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(log.weakdep.res,file="log_weakdep_quantiles.Rdata")}
  
  # Gaussian
  set.seed(k)
  x<-qexp(pnorm(rmvnorm(5000,sigma=matrix(c(1,0.8,0.8,1),2,2))))
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  # qr = QR.2d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  gauss.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(gauss.res,file="gauss_quantiles.Rdata")}
  
  # inverted logistic
  set.seed(k)
  x<-rbvevd(5000,dep=0.7,mar1=c(1,1,1))
  x<-1/x
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  tau=0.95
  # qr = QR.2d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  invlog.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(invlog.res,file="invlog_quantiles.Rdata")}
}
