rm(list=ls())

library(mvtnorm)
library(evd)
library(qgam)

# Conditional extremes calculations:
# setwd("~/SimStudy2d_code/")
setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_2d/simstudypar/")
source("Extra Functions/CEproflik.r")

library(geometricMVE)
#======================
# To install for the first time from the current working directory:
# library(devtools)
# install_local("geometricMVE_0.1.5.tar.gz")
#======================

# file.sources = list.files("~/Dropbox/phd_research/geomMVE_work/r_code/geometricMVE_0.1.5/R", 
#                           pattern="*.R$", full.names=TRUE, 
#                           ignore.case=TRUE)
# sapply(file.sources,source,.GlobalEnv)
#=================================

# Read in spectral densities and NLL functions
source("Extra Functions/spectral_densities.R")
spec.dens.fns = list(h.log,h.neglog,h.bilog,h.negbilog,h.ct)
H.fns = list(H.log.2d,H.neglog.2d,H.bilog.2d,H.negbilog.2d,H.ct.2d)
nll.fns = c(nll.log,nll.neglog,nll.bilog,nll.negbilog,nll.ct)
starting.vals = list(c(0.5),c(0.5),rep(0.5,2),rep(0.5,2),rep(0.5,2))


#=================================

nrep<-200
llvals<-matrix(0,ncol=4,nrow=nrep)
parsA<-matrix(0,ncol=2,nrow=nrep)
parsB<-matrix(0,ncol=2,nrow=nrep)
parsC<-matrix(0,ncol=2,nrow=nrep)
parsD<-matrix(0,ncol=2,nrow=nrep)

llvals2<-matrix(0,ncol=4,nrow=nrep)
pars2A<-matrix(0,ncol=2,nrow=nrep)
pars2B<-matrix(0,ncol=2,nrow=nrep)
pars2C<-matrix(0,ncol=2,nrow=nrep)
pars2D<-matrix(0,ncol=2,nrow=nrep)

etahat<-NULL
empgauge<-empgauge2<-list()
kvalsS<-kvalsE<-matrix(0,nrow=nrep,ncol=2)

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

v1<-7

x110<-x11-v1
x120<-x12-v1
y110<-y11-v1
y120<-y12-v1

# Set 2------------------

x21<-10
x22<-12
y21<-6
y22<-8

v2<-6

x210<-x21-v2
x220<-x22-v2
y210<-y21-v2
y220<-y22-v2

# Set 3------------------

x31<-10
x32<-12
y31<-2
y32<-4

v3<-2

x310<-x31-v3
x320<-x32-v3
y310<-y31-v3
y320<-y32-v3


prob1.est.geomE<-prob1.est.geomE2<-prob1.est.geomS<-prob1.est.geomS2<-prob1.est.MRV<-prob1.est.HRV<-prob1.est.CE<-prob1.est.CE2<-numeric(nrep)
prob2.est.geomE<-prob2.est.geomE2<-prob2.est.geomS<-prob2.est.geomS2<-prob2.est.MRV<-prob2.est.HRV<-prob2.est.CE<-prob2.est.CE2<-numeric(nrep)
prob3.est.geomE<-prob3.est.geomE2<-prob3.est.geomS<-prob3.est.geomS2<-prob3.est.MRV<-prob3.est.HRV<-prob3.est.CE<-prob3.est.CE2<-numeric(nrep)

prob1.est.geomE3<-prob1.est.geomS3<-NULL

prob1.est.parMRV.samp<-prob2.est.parMRV.samp<-prob3.est.parMRV.samp<-prob1.est.parMRV.int<-prob2.est.parMRV.int<-prob3.est.parMRV.int<-numeric(nrep)

aAI<-0.7
n<-5000
n1<-50000

prob1.true = true.prob.in.B.d2(x11,x12,y11,y12,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)
prob2.true = true.prob.in.B.d2(x21,x22,y21,y22,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)
prob3.true = true.prob.in.B.d2(x31,x32,y31,y32,
                               distn.fn=inv.logistic.cdf,dep.par=aAI)

time.start<-Sys.time()
# set.seed(321)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  set.seed(k)
  x<-rbvevd(n,dep=aAI,mar1=c(1,1,1))
  x<-1/x  
   
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  
  #==================================
  # 1. Using qgam
  
  qr<-QR.2d(r=r,w=w,method = "smooth")
  
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind]
  r0w<-qr$r0w[excind]
  
  ### For empirical gauge estimation
  
  empgauge[[k]]<-cbind(qr$r.tau.wpts*qr$wpts/max(qr$r.tau.wpts*qr$wpts),qr$r.tau.wpts*(1-qr$wpts)/max(qr$r.tau.wpts*(1-qr$wpts)))
  
  ##################################
  
  fitA<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,gfun = gauge_rvad,init.val = c(2,0.5), 
                            pos.par=TRUE, lower.limit=c(0,0),upper.limit=c(1000,1))
  fitB<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,gfun = gauge_gaussian,init.val = c(2,0.5),
                            pos.par=TRUE, lower.limit=c(0,-1),upper.limit=c(1000,1))
  fitC<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,gfun = gauge_invlogistic,init.val = c(2,2),
                            pos.par=TRUE, lower.limit=c(0,1),upper.limit=c(1000,5000))
  fitD<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,gfun = gauge_square,init.val = c(2,0.5),
                            pos.par=TRUE, lower.limit=c(0,0),upper.limit=c(1000,1))
  
  parsA[k,]<-fitA$mle
  parsB[k,]<-fitB$mle
  parsC[k,]<-fitC$mle
  parsD[k,]<-fitD$mle
  
  llvals[k,]<-c(fitA$nllh,fitB$nllh,fitC$nllh,fitD$nllh)
  
  if(which.min(llvals[k,])==1){fittedpar<-fitA$mle;fittedgfun<-gauge_rvad}
  if(which.min(llvals[k,])==2){fittedpar<-fitB$mle;fittedgfun<-gauge_gaussian}
  if(which.min(llvals[k,])==3){fittedpar<-fitC$mle;fittedgfun<-gauge_invlogistic}
  if(which.min(llvals[k,])==4){fittedpar<-fitD$mle;fittedgfun<-gauge_square}
  
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fittedpar,gfun=fittedgfun)
  
  prob1.est.geomS[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)
  prob2.est.geomS[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * length(rexc)/length(r)
  prob3.est.geomS[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * length(rexc)/length(r)
  
  ## Extrapolation from higher levels
  
  # Find appropriate values of k (for R'>k)
  
  cm<-max(qr$wpts*qr$r.tau.wpts,max((1-qr$wpts)*qr$r.tau.wpts))
  km<-10/cm #  max value that works everywhere given the regions we're interested in
  
  
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fittedpar,gfun=fittedgfun,k=km)
  
  prob1.est.geomS2[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob2.est.geomS2[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob3.est.geomS2[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * Rexc.prob.k.2d(k=km, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  
  # Higher level just for set 1
  
  cs<-NULL
  kgrid<-seq(km,3.8,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.2d(x=c(x11,x12),y=c(y11,y12),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  
  ks<-max(kgrid[cs])
  
  xstar<-sim.2d(w=wexc,r0w = r0w,nsim = n1,par=fittedpar,gfun=fittedgfun,k=ks)
  prob1.est.geomS3[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * Rexc.prob.k.2d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  
  kvalsS[k,]<-c(km,ks)
  
  #==================================
  # 2. Using empirical estimation
  
  ########################
  
  qr<-QR.2d(r=r,w=w,method = "empirical")
  
  excind<-r>qr$r0w
  rexc2<-r[excind]
  wexc2<-w[excind]
  r0w2<-qr$r0w[excind]
  
  ### For empirical gauge estimation
  
  empgauge2[[k]]<-cbind(qr$r.tau.wpts*qr$wpts/max(qr$r.tau.wpts*qr$wpts),qr$r.tau.wpts*(1-qr$wpts)/max(qr$r.tau.wpts*(1-qr$wpts)))
  
  if(any(is.na(r0w2))){
    na.ind<-which(is.na(r0w2))
    rexc2<-rexc2[-na.ind]
    wexc2<-wexc2[-na.ind]
    r0w2<-r0w2[-na.ind]
  }
  
  fit2A<-fit.geometricMVE.2d(r=rexc2,w=wexc2,r0w=r0w2,gfun = gauge_rvad,init.val = c(2,0.5),
                             pos.par=TRUE, lower.limit=c(0,0),upper.limit=c(1000,1))
  fit2B<-fit.geometricMVE.2d(r=rexc2,w=wexc2,r0w=r0w2,gfun = gauge_gaussian,init.val = c(2,0.5),
                             pos.par=TRUE, lower.limit=c(0,-1),upper.limit=c(1000,1))
  fit2C<-fit.geometricMVE.2d(r=rexc2,w=wexc2,r0w=r0w2,gfun = gauge_invlogistic,init.val = c(2,2),
                             pos.par=TRUE, lower.limit=c(0,1),upper.limit=c(1000,5000))
  fit2D<-fit.geometricMVE.2d(r=rexc2,w=wexc2,r0w=r0w2,gfun = gauge_square,init.val = c(2,0.5),
                             pos.par=TRUE, lower.limit=c(0,0),upper.limit=c(1000,1))
  
  pars2A[k,]<-fit2A$mle
  pars2B[k,]<-fit2B$mle
  pars2C[k,]<-fit2C$mle
  pars2D[k,]<-fit2D$mle
  
  llvals2[k,]<-c(fit2A$nllh,fit2B$nllh,fit2C$nllh,fit2D$nllh)
  
  if(which.min(llvals2[k,])==1){fittedpar<-fit2A$mle;fittedgfun<-gauge_rvad}
  if(which.min(llvals2[k,])==2){fittedpar<-fit2B$mle;fittedgfun<-gauge_gaussian}
  if(which.min(llvals2[k,])==3){fittedpar<-fit2C$mle;fittedgfun<-gauge_invlogistic}
  if(which.min(llvals2[k,])==4){fittedpar<-fit2D$mle;fittedgfun<-gauge_square}
  
  xstar2<-sim.2d(w=wexc2,r0w = r0w2,nsim = n1,par=fittedpar,gfun=fittedgfun)
  
  prob1.est.geomE[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12) * length(rexc2)/length(r)
  prob2.est.geomE[k]<-mean(xstar2[,1]>x21&xstar2[,1]<x22&xstar2[,2]>y21&xstar2[,2]<y22) * length(rexc2)/length(r)
  prob3.est.geomE[k]<-mean(xstar2[,1]>x31&xstar2[,1]<x32&xstar2[,2]>y31&xstar2[,2]<y32) * length(rexc2)/length(r)
  
  
  # Higher level
  cm<-max(qr$wpts*qr$r.tau.wpts,max((1-qr$wpts)*qr$r.tau.wpts))
  km<-10/cm #  max value that works everywhere given the regions we're interested in
  
  xstar2<-sim.2d(w=wexc2,r0w = r0w2,nsim = n1,par=fittedpar,gfun=fittedgfun,k=km)
  
  prob1.est.geomE2[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12) * Rexc.prob.k.2d(k=km, r0w=r0w2, w=wexc2, gfun=fittedgfun, par=fittedpar)* length(rexc2)/length(r)
  prob2.est.geomE2[k]<-mean(xstar2[,1]>x21&xstar2[,1]<x22&xstar2[,2]>y21&xstar2[,2]<y22) * Rexc.prob.k.2d(k=km, r0w=r0w2, w=wexc2, gfun=fittedgfun, par=fittedpar)* length(rexc2)/length(r)
  prob3.est.geomE2[k]<-mean(xstar2[,1]>x31&xstar2[,1]<x32&xstar2[,2]>y31&xstar2[,2]<y32) * Rexc.prob.k.2d(k=km, r0w=r0w2, w=wexc2, gfun=fittedgfun, par=fittedpar)* length(rexc2)/length(r)
  
  # Higher level just for set 3
  cs<-NULL
  kgrid<-seq(km,3.8,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.2d(x=c(x11,x12),y=c(y11,y12),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  
  ks<-max(kgrid[cs])
  

  xstar2<-sim.2d(w=wexc2,r0w = r0w2,nsim = n1,par=fittedpar,gfun=fittedgfun,k=ks)
  prob1.est.geomE3[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12) * Rexc.prob.k.2d(k=ks, r0w=r0w2, w=wexc2, gfun=fittedgfun, par=fittedpar)* length(rexc2)/length(r)
 
  kvalsE[k,]<-c(km,ks)
  
  #==================================
  
  ########################
  
  # HT comparison 
  th<-quantile(x[,1],0.95)
  f1<-optim(HT.nll.PL,par=c(0.5,0.5),dat=x,u=th,which=1)
  Z1<-(x[x[,1]>th,2]-f1$par[1]*x[x[,1]>th,1])/x[x[,1]>th,1]^f1$par[2]
  xe.star<-rexp(n1)+qexp(0.999)
  Z1star<-sample(Z1,size=n1,replace=T)
  Y1.star<-f1$par[1]*xe.star+Z1star*xe.star^f1$par[2]
  xstar<-cbind(xe.star,Y1.star)
  prob1.est.CE[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * 0.001
  prob2.est.CE[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * 0.001
  prob3.est.CE[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * 0.001
  
  # Higher level
  xe.star<-rexp(n1)+10
  Z1star<-sample(Z1,size=n1,replace=T)
  Y1.star<-f1$par[1]*xe.star+Z1star*xe.star^f1$par[2]
  xstar<-cbind(xe.star,Y1.star)
  prob1.est.CE2[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * exp(-10)
  prob2.est.CE2[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * exp(-10)
  prob3.est.CE2[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * exp(-10)
  
  
  #HRV
  xmin<-apply(x,1,min)
  etahat[k]<-mean(xmin[xmin>quantile(xmin,0.95)]-quantile(xmin,0.95))
  
  
  prob1.est.HRV[k]<-mean(x[,1]>x110&x[,1]<x120&x[,2]>y110&x[,2]<y120)*exp(-v1/etahat[k])
  prob2.est.HRV[k]<-mean(x[,1]>x210&x[,1]<x220&x[,2]>y210&x[,2]<y220)*exp(-v2/etahat[k])
  prob3.est.HRV[k]<-mean(x[,1]>x310&x[,1]<x320&x[,2]>y310&x[,2]<y320)*exp(-v3/etahat[k])
  
  # MRV
  prob1.est.MRV[k]<-mean(x[,1]>x110&x[,1]<x120&x[,2]>y110&x[,2]<y120)*exp(-v1)
  prob2.est.MRV[k]<-mean(x[,1]>x210&x[,1]<x220&x[,2]>y210&x[,2]<y220)*exp(-v2)
  prob3.est.MRV[k]<-mean(x[,1]>x310&x[,1]<x320&x[,2]>y310&x[,2]<y320)*exp(-v3)
  
  # Parametric MRV
  x.P = exp(x)
  
  r<-apply(x.P,1,sum)
  w<-x.P[,1]/r
  r.thresh = quantile(r,p=0.95,names=F)
  
  excind = r>r.thresh
  rexc = r[excind]
  wexc = w[excind]
  
  optim.res = mapply(function(f,start){
    res = try(optim(par=start,fn=f,wexc=wexc))
    while(class(res)=="try-error"){
      # res = list(par=start,value=Inf)
      alt.start = runif(length(start))
      res = try(optim(par=alt.start,fn=f,wexc=wexc))
    }
    res$aic = 2*(length(res$par) + res$value)
    return(res)
  },
  f=nll.fns,start=starting.vals)
  # if(all(class(optim.res)=="list")){
  #   optim.res=do.call(cbind,optim.res)
  # }
  # optim.res = mapply(function(f,start){
  #   res = optim(par=start,fn=f,wexc=wexc)
  #   res$aic = 2*(length(res$par) + res$value)
  #   return(res)
  # },
  # f=nll.fns,start=starting.vals)
  
  best.h.idx = which.min(optim.res["aic",])
  best.h = spec.dens.fns[[best.h.idx]]
  best.H = H.fns[[best.h.idx]]
  h.mle = optim.res["par",][[best.h.idx]]
  
  ## version 1 - sampling-based approach
  w.mesh = seq(0,1,length.out=1000)
  M = max((best.h(w = w.mesh, dep = h.mle))/
            dbeta(w.mesh,0.5,0.5), na.rm=T) + 1
  wstar = numeric(n1)
  idx = 1
  iter.count = 1
  max.iters = 1e6
  max.reached = F
  while(idx<=n1){
    u = runif(1)
    prop =  rbeta(1,0.5,0.5)
    acc.rej.cond = (best.h(w = prop, dep = h.mle)) / (M*dbeta(prop,0.5,0.5))
    if(u < acc.rej.cond){
      wstar[idx] = prop
      idx = idx + 1
    }
    iter.count = iter.count + 1
    if(iter.count > max.iters){
      message("Max iters reached in rejection sampling")
      prob1.est.parMRV.samp[k] = NA
      prob2.est.parMRV.samp[k] = NA
      prob3.est.parMRV.samp[k] = NA
      max.reached = T
    }
  }
  
  if(!max.reached){
    rstar = 1/(1-runif(n1))
    c1 = sum(exp(c(x11,y11)))#9*n
    xstar1.P = (rstar*c1)*cbind(wstar,1-wstar)
    c2 = sum(exp(c(x21,y21)))#4*n
    xstar2.P = (rstar*c2)*cbind(wstar,1-wstar)
    c3 = sum(exp(c(x31,y31)))#4*n
    xstar3.P = (rstar*c3)*cbind(wstar,1-wstar)
    
    prob1.est.parMRV.samp[k] = mean(excind)*(r.thresh/c1)*#(1/c1)*
      mean(xstar1.P[,1] > exp(x11) & xstar1.P[,1] < exp(x12) &
             xstar1.P[,2] > exp(y11) & xstar1.P[,2] < exp(y12))
    prob2.est.parMRV.samp[k] = mean(excind)*(r.thresh/c2)*#(1/c2)*
      mean(xstar2.P[,1] > exp(x21) & xstar2.P[,1] < exp(x22) &
             xstar2.P[,2] > exp(y21) & xstar2.P[,2] < exp(y22))
    prob3.est.parMRV.samp[k] = mean(excind)*(r.thresh/c3)*
      mean(xstar3.P[,1] > exp(x31) & xstar3.P[,1] < exp(x32) &
             xstar3.P[,2] > exp(y31) & xstar3.P[,2] < exp(y32))
  }
  
  ## version 2 - integration-based approach
  integrand = function(s,c,dep.pars,xl,xu,yl,yu){
    pmax((best.H(min(xu/s,1-(yl/s)),dep=dep.pars) - 
            best.H(max(xl/s,1-(yu/s)),dep=dep.pars)),0)*
      (c/(s^2))
  }
  integrand.vec = Vectorize(integrand,vectorize.args = "s")
  
  prob1.est.parMRV.int[k] = mean(excind)*(r.thresh/c1)*
    integrate(f=integrand.vec,
              lower=c1,
              upper=sum((exp(c(x12,y12)))),
              c=c1,
              dep.pars=h.mle,
              xl=exp(x11),xu=exp(x12),yl=exp(y11),yu=exp(y12))$value
  prob2.est.parMRV.int[k] = mean(excind)*(r.thresh/c2)*
    integrate(f=integrand.vec,
              lower=c2,
              upper=sum((exp(c(x22,y22)))),
              c=c2,
              dep.pars=h.mle,
              xl=exp(x21),xu=exp(x22),yl=exp(y21),yu=exp(y22))$value
  prob3.est.parMRV.int[k] = mean(excind)*(r.thresh/c3)*
    integrate(f=integrand.vec,
              lower=c3,
              upper=sum((exp(c(x32,y32)))),
              c=c3,
              dep.pars=h.mle,
              xl=exp(x31),xu=exp(x32),yl=exp(y31),yu=exp(y32))$value
  
  print(k)
  if(k%%20==0){save.image("SimStudy2d_invlog.Rdata")}
}

time.end<-Sys.time()
save.image("SimStudy2d_invlog.Rdata")
#save.image("~/Dropbox/phd_research/geomMVE_work/r_code/SimStudy2d_code/SimStudy2d_invlog.Rdata")


#load("~/Dropbox/phd_research/geomMVE_work/r_code/SimStudy2d_code/SimStudy2d_invlog.Rdata")
