rm(list=ls())

library(mvtnorm)
library(evd)
library(qgam)

# fixed shape parameter in truncated gamma
source("~/Dropbox/phd_research/geomMVE_work/r_code/Functions/fit_geometricMVE_fixedshape.R")
source("~/Dropbox/phd_research/geomMVE_work/r_code/Functions/rw_gamma_lik_fixedshape.R")

setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_3d/simstudypar2/")

library(geometricMVE)

#===========================================

# Define extreme sets of interest, and the true probabilities of lying within them

x11<-8
x12<-10
y11<-8
y12<-10

x21<-8
x22<-10
y21<-5
y22<-7

x31<-8
x32<-10
y31<-2
y32<-4

z11<-0.01
z12<-3
z21<-0.01
z22<-3
z31<-0.01
z32<-3

asy <- list(.0, .0, .0, c(.5,0.5), c(.5,0.5), c(.5,.5), c(.0,.0,.0))

p1<-pmvevd(c(qfrechet(pexp(c(x12,y12,z12)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p2<-pmvevd(c(qfrechet(pexp(c(x11,y12,z12)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p3<-pmvevd(c(qfrechet(pexp(c(x12,y11,z12)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p4<-pmvevd(c(qfrechet(pexp(c(x12,y12,z11)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p5<-pmvevd(c(qfrechet(pexp(c(x11,y11,z12)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p6<-pmvevd(c(qfrechet(pexp(c(x11,y12,z11)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p7<-pmvevd(c(qfrechet(pexp(c(x12,y11,z11)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p8<-pmvevd(c(qfrechet(pexp(c(x11,y11,z11)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))

TP1<-p1-p2-p3-p4+p5+p6+p7-p8

p1<-pmvevd(c(qfrechet(pexp(c(x22,y22,z22)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p2<-pmvevd(c(qfrechet(pexp(c(x21,y22,z22)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p3<-pmvevd(c(qfrechet(pexp(c(x22,y21,z22)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p4<-pmvevd(c(qfrechet(pexp(c(x22,y22,z21)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p5<-pmvevd(c(qfrechet(pexp(c(x21,y21,z22)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p6<-pmvevd(c(qfrechet(pexp(c(x21,y22,z21)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p7<-pmvevd(c(qfrechet(pexp(c(x22,y21,z21)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p8<-pmvevd(c(qfrechet(pexp(c(x21,y21,z21)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))

TP2<-p1-p2-p3-p4+p5+p6+p7-p8

p1<-pmvevd(c(qfrechet(pexp(c(x32,y32,z32)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p2<-pmvevd(c(qfrechet(pexp(c(x31,y32,z32)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p3<-pmvevd(c(qfrechet(pexp(c(x32,y31,z32)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p4<-pmvevd(c(qfrechet(pexp(c(x32,y32,z31)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p5<-pmvevd(c(qfrechet(pexp(c(x31,y31,z32)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p6<-pmvevd(c(qfrechet(pexp(c(x31,y32,z31)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p7<-pmvevd(c(qfrechet(pexp(c(x32,y31,z31)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
p8<-pmvevd(c(qfrechet(pexp(c(x31,y31,z31)),0,1,1)), dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))

TP3<-p1-p2-p3-p4+p5+p6+p7-p8

#==========================================
n1<-50000
nrep<-200

llvals<-AICvals<-matrix(NA,ncol=3,nrow=nrep)
prob1.est.geom<-prob2.est.geom<-prob3.est.geom<-NULL
prob1.est.geom2<-prob2.est.geom2<-prob3.est.geom2<-NULL
prob1.est.geom3<-prob2.est.geom3<-prob3.est.geom3<-NULL

Set1inRegion1<-Set2inRegion1<-Set3inRegion1<-numeric(nrep)

load("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_3d/quantiles/alog1_quantiles.Rdata")

time.start<-Sys.time()
# set.seed(101)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  set.seed(k)
  
  x<-alog1.res[[k]]$x
  r<-apply(x,1,sum)
  w<-x/r
  
  qr<-alog1.res[[k]]$qr
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind,]
  r0w<-qr$r0w[excind]
  na.ind<-which(is.na(rexc))
  na.ind
  if(length(na.ind)>0){
    rexc<-rexc[-na.ind]
    wexc<-wexc[-na.ind,]
    r0w<-r0w[-na.ind]}
  
  # Check if sets are in regions R'>1, R'>2
  Set1inRegion1[k]<-checkset.3d(x=c(x11,x12),y=c(y11,y12),z=c(z11,z12),w=w,r0w=qr$r0w)
  Set2inRegion1[k]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),w=w,r0w=qr$r0w)
  Set3inRegion1[k]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),w=w,r0w=qr$r0w)

  fit<-fit.geometricMVE.3d.fixedshape(r = rexc, w=wexc, r0w = r0w, gfun = gauge_rvad_pw3d,
                                      init.val = rep(0.5,3),
                                      pos.par=T,fixed.shape=T,
                                      control=list(maxit=1e5))
  
  #print(fit)
  #stop()
  
  # llvals[k,]<-fit$nllh

  fittedpar=fit$mle
  fittedgfun=gauge_rvad_pw3d
  
  xstar<-sim.3d(w = wexc,r0w = r0w,nsim = n1,par = fittedpar, gfun = fittedgfun)
  
  prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12&xstar[,3]>z11&xstar[,3]<z12) * length(rexc)/length(r)
  prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22&xstar[,3]>z21&xstar[,3]<z22) * length(rexc)/length(r)
  prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32&xstar[,3]>z31&xstar[,3]<z32) * length(rexc)/length(r)
  
  # Simulate above higher thresholds (for sets 2 and 3)
  cm<-max(max(c(qr$wpts[,1]*qr$r.tau.wpts),na.rm=T),max(qr$wpts[,2]*qr$r.tau.wpts,na.rm=T),max((1-qr$wpts[,1]-qr$wpts[,2])*qr$r.tau.wpts, na.rm=T))
  km<-8/cm #  max value that works everywhere given the regions we're interested in
  
  # Higher level for set 2
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  ks=ifelse(ks<1,1,ks)
  xstar2<-sim.3d(k=ks,w = wexc,r0w = r0w,nsim = n1,par = fittedpar, gfun = fittedgfun)
 
  prob1.est.geom2[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12&xstar2[,3]>z11&xstar2[,3]<z12) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob2.est.geom2[k]<-mean(xstar2[,1]>x21&xstar2[,1]<x22&xstar2[,2]>y21&xstar2[,2]<y22&xstar2[,3]>z21&xstar2[,3]<z22) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob3.est.geom2[k]<-mean(xstar2[,1]>x31&xstar2[,1]<x32&xstar2[,2]>y31&xstar2[,2]<y32&xstar2[,3]>z31&xstar2[,3]<z32) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  
  # Higher level just for set 3
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  ks=ifelse(ks<1,1,ks)
  
  xstar3<-sim.3d(k=ks,w = wexc,r0w = r0w,nsim = n1,par = fittedpar, gfun = fittedgfun)
 
  prob1.est.geom3[k]<-mean(xstar3[,1]>x11&xstar3[,1]<x12&xstar3[,2]>y11&xstar3[,2]<y12&xstar3[,3]>z11&xstar3[,3]<z12) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob2.est.geom3[k]<-mean(xstar3[,1]>x21&xstar3[,1]<x22&xstar3[,2]>y21&xstar3[,2]<y22&xstar3[,3]>z21&xstar3[,3]<z22) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob3.est.geom3[k]<-mean(xstar3[,1]>x31&xstar3[,1]<x32&xstar3[,2]>y31&xstar3[,2]<y32&xstar3[,3]>z31&xstar3[,3]<z32) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  
  ###################################
  if(k%%20==0){save.image("SimStudy3d_alog1.Rdata")}
}
rm(alog1.res)
time.end<-Sys.time()

save.image("SimStudy3d_alog1.Rdata")
