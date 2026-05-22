rm(list=ls())

library(mvtnorm)
library(evd)
library(qgam)

# fixed shape parameter in truncated gamma
source("~/Dropbox/phd_research/geomMVE_work/r_code/Functions/fit_geometricMVE_fixedshape.R")
source("~/Dropbox/phd_research/geomMVE_work/r_code/Functions/rw_gamma_lik_2d_fixedshape.R")

# Conditional extremes calculations:
setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/par/simstudypar2/")

library(geometricMVE)

#----------------------------------

nrep<-200
llvals<-numeric(nrep)
pars<-matrix(0,ncol=2,nrow=nrep)

inv.logistic.cdf <- function(x.vec,dep.par=0.5){
  # input - vector in Exp(1) margins
  1 - sum(exp(-x.vec)) + exp(-(sum(x.vec^(1/dep.par)))^dep.par)
}
true.prob.in.B.d2 <- function(xmin,xmax,ymin,ymax,distn.fn,...){
  # inputs in Exp(1) margins
  return(distn.fn(c(xmax,ymax),...) - distn.fn(c(xmax,ymin),...) - distn.fn(c(xmin,ymax),...) +
           distn.fn(c(xmin,ymin),...))
}

# Set 1-----------------

x11<-10
x12<-12
y11<-10
y12<-12

# Set 2------------------

x21<-10
x22<-12
y21<-6
y22<-8

# Set 3------------------

x31<-10
x32<-12
y31<-2
y32<-4

prob1.est.geom<-prob1.est.geom2<-numeric(nrep)
prob2.est.geom<-prob2.est.geom2<-numeric(nrep)
prob3.est.geom<-prob3.est.geom2<-numeric(nrep)
prob1.est.geom3<-NULL

aAI<-0.7
n<-5000
n1<-50000

prob1.true = true.prob.in.B.d2(x11,x12,y11,y12,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)
prob2.true = true.prob.in.B.d2(x21,x22,y21,y22,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)
prob3.true = true.prob.in.B.d2(x31,x32,y31,y32,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)

load("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/quantiles/quantiles_kde/invlog_quantiles.Rdata")

time.start<-Sys.time()
# set.seed(321)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  set.seed(k)
  x<-invlog.res[[k]]$x
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  
  #==================================
  
  qr<-invlog.res[[k]]$qr
  
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind]
  r0w<-qr$r0w[excind]
  
  fit = fit.geometricMVE.2d.fixedshape(r=rexc,w=wexc,r0w=r0w,gfun=gauge_invlogistic,
                                       init.val=0.5,fixed.shape=T,pos.par=T)
  pars[k,]<-fit$mle
  llvals[k]<-fit$nllh
  
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fit$mle,gfun=gauge_invlogistic)
  
  prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)
  prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * length(rexc)/length(r)
  prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * length(rexc)/length(r)
  
  ## Extrapolation from higher levels
  # Find appropriate values of k (for R'>k)
  cm<-max(max((qr$wpts)*qr$r.tau.wpts),max((1-qr$wpts)*qr$r.tau.wpts))
  km<-10/cm #  max value that works everywhere given the regions we're interested in
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fit$mle,gfun=gauge_invlogistic,k=km)
  
  prob1.est.geom2[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=gauge_invlogistic, par=fit$mle)* length(rexc)/length(r)
  prob2.est.geom2[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=gauge_invlogistic, par=fit$mle)* length(rexc)/length(r)
  prob3.est.geom2[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=gauge_invlogistic, par=fit$mle)* length(rexc)/length(r)
  
  # Higher level just for set 1
  cs<-NULL
  kgrid<-seq(km,4.2,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.2d(x=c(x11,x12),y=c(y11,y12),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fit$mle,gfun=gauge_invlogistic,k=ks)
  prob1.est.geom3[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * Rexc.prob.k.2d(k=ks, r0w=r0w, w=wexc, gfun=gauge_invlogistic, par=fit$mle)* length(rexc)/length(r)
  
  if(k%%20==0){save.image("SimStudy2d_invlog.Rdata")}
}

rm(invlog.res)
time.end<-Sys.time()
save.image("SimStudy2d_invlog.Rdata")
