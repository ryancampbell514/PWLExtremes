rm(list=ls(all=T))

library(evd)
library(geometricMVE)
library(interp)
library(mvtnorm)

setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_3d/quantiles_emp/")
# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

alog1.res = NULL
alog2.res = NULL
mix.res = NULL

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
  qr = geometricMVE::QR.3d(r,w,tau=tau)
  alog1.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(alog1.res,file="alog1_quantiles.Rdata")}
  
  # alog2
  set.seed(k)
  asy <- list(0.5, 0, 0, c(0.5,0.5), c(0,0), c(0.5,1), c(.0,.0,.0))
  x<-rmvevd(5000, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
  x<--log(1-exp(-1/x))
  r<-apply(x,1,sum)
  w<-x/r
  tau=0.95
  qr = geometricMVE::QR.3d(r,w,tau=tau)
  alog2.res[k] = mean(r>qr$r0w,na.rm=T)
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
  qr = geometricMVE::QR.3d(r,w,tau=tau)
  mix.res[k] = mean(r>qr$r0w,na.rm=T)
  if(k%%20==0){save(mix.res,file="mix_quantiles.Rdata")}
}
