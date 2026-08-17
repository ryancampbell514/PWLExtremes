rm(list=ls())

library(evd)
# library(qgam)
library(rmutil)
library(geometry)
library(scales)
# library(geometricMVE)
library(mvtnorm)
library(interp)

# par(pty="s")
# setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/quantiles/")
# 
# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# library(this.path)
# setwd(this.path::here())
setwd("/home/rcampbell/lancaster/pw_lin_gauge/SimStudy_3d/quantiles/quantiles_kde")

fn.dir = "../../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# source("../../../kde_quants_bw_plots_unconstrained/thresh_estimation_unconstrained.R")
source("../../../offsimplex_threshold/thresh_estimation_constrained_AKDE.R")

alog1.res = list(list(NULL))
alog2.res = list(list(NULL))
mix.res = list(list(NULL))

nrep=200
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  # alog1
  set.seed(k)
  asy <- list(.0, .0, .0, c(.5,0.5), c(.5,0.5), c(.5,.5), c(.0,.0,.0))
  x<-rmvevd(5000, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
  x<--log(1-exp(-1/x))
  r<-apply(x,1,sum)
  w<-x/r
  tau=0.95
  # qr = QR.3d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  alog1.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(alog1.res,file="alog1_quantiles.Rdata")}
  
  # alog2
  set.seed(k)
  asy <- list(0.5, 0, 0, c(0.5,0.5), c(0,0), c(0.5,1), c(.0,.0,.0))
  x<-rmvevd(5000, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
  x<--log(1-exp(-1/x))
  r<-apply(x,1,sum)
  w<-x/r
  tau=0.95
  # qr = QR.3d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  alog2.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(alog2.res,file="alog2_quantiles.Rdata")}
  
  # mixture model
  set.seed(k)
  asy <- list(0, 0, 0, c(0.5,0.5), c(0,0), c(0,0), c(0.5,0.5,1))
  rho=0.4
  phi=0.6
  x.norm = qexp(pnorm(rmvnorm(2500,sigma=matrix(c(1,rho,rho,
                                                  rho,1,rho,
                                                  rho,rho,1),3,3))))
  x.alog = rmvevd(2500, dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1))
  x.alog<--log(1-exp(-1/x.alog))
  x = rbind(x.norm,x.alog)
  r<-apply(x,1,sum)
  w<-x/r
  tau=0.95
  # qr = QR.3d.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)
  qr = fit.thresh(r,w,tau=tau,method="KDE")
  mix.res[[k]] = list(x=x,qr=qr)
  if(k%%20==0){save(mix.res,file="mix_quantiles.Rdata")}
}
