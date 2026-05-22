rm(list=ls(all=T))

library(evd)
library(mvtnorm)

# Conditional extremes and tau calculations:
setwd("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_3d/simstudypar/")
source("Extra Functions/CEproflik.r")
source("Extra Functions/SimpsonEtAlFunctions.R")

library(geometricMVE)
#======================
# To install for the first time from the current working directory:
# library(devtools)
# install_local("geometricMVE_0.1.5.tar.gz")
#======================


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

asy <- list(0, 0, 0, c(0.5,0.5), c(0,0), c(0,0), c(0.5,0.5,1))
rho=0.4
phi=0.6

p1<-0.5*pmvevd(c(qfrechet(pexp(c(x12,y12,z12)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x12,y12,z12))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p2<-0.5*pmvevd(c(qfrechet(pexp(c(x11,y12,z12)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x11,y12,z12))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p3<-0.5*pmvevd(c(qfrechet(pexp(c(x12,y11,z12)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x12,y11,z12))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p4<-0.5*pmvevd(c(qfrechet(pexp(c(x12,y12,z11)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x12,y12,z11))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p5<-0.5*pmvevd(c(qfrechet(pexp(c(x11,y11,z12)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x11,y11,z12))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p6<-0.5*pmvevd(c(qfrechet(pexp(c(x11,y12,z11)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x11,y12,z11))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p7<-0.5*pmvevd(c(qfrechet(pexp(c(x12,y11,z11)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x12,y11,z11))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p8<-0.5*pmvevd(c(qfrechet(pexp(c(x11,y11,z11)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x11,y11,z11))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]

TP1<-p1-p2-p3-p4+p5+p6+p7-p8

p1<-0.5*pmvevd(c(qfrechet(pexp(c(x22,y22,z22)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x22,y22,z22))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p2<-0.5*pmvevd(c(qfrechet(pexp(c(x21,y22,z22)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x21,y22,z22))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p3<-0.5*pmvevd(c(qfrechet(pexp(c(x22,y21,z22)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x22,y21,z22))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p4<-0.5*pmvevd(c(qfrechet(pexp(c(x22,y22,z21)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x22,y22,z21))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p5<-0.5*pmvevd(c(qfrechet(pexp(c(x21,y21,z22)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x21,y21,z22))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p6<-0.5*pmvevd(c(qfrechet(pexp(c(x21,y22,z21)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x21,y22,z21))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p7<-0.5*pmvevd(c(qfrechet(pexp(c(x22,y21,z21)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x22,y21,z21))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p8<-0.5*pmvevd(c(qfrechet(pexp(c(x21,y21,z21)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x21,y21,z21))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]

TP2<-p1-p2-p3-p4+p5+p6+p7-p8

p1<-0.5*pmvevd(c(qfrechet(pexp(c(x32,y32,z32)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x32,y32,z32))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p2<-0.5*pmvevd(c(qfrechet(pexp(c(x31,y32,z32)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x31,y32,z32))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p3<-0.5*pmvevd(c(qfrechet(pexp(c(x32,y31,z32)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x32,y31,z32))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p4<-0.5*pmvevd(c(qfrechet(pexp(c(x32,y32,z31)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x32,y32,z31))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p5<-0.5*pmvevd(c(qfrechet(pexp(c(x31,y31,z32)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x31,y31,z32))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p6<-0.5*pmvevd(c(qfrechet(pexp(c(x31,y32,z31)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x31,y32,z31))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p7<-0.5*pmvevd(c(qfrechet(pexp(c(x32,y31,z31)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x32,y31,z31))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]
p8<-0.5*pmvevd(c(qfrechet(pexp(c(x31,y31,z31)),0,1,1)), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1)) +
  0.5*pmvnorm(upper=qnorm(pexp(c(x31,y31,z31))),sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),3,3))[1]

TP3<-p1-p2-p3-p4+p5+p6+p7-p8

#==========================================
n1<-50000
nrep<-200

llvals<-AICvals<-matrix(NA,ncol=4,nrow=nrep)
prob1.est.geom<-prob2.est.geom<-prob3.est.geom<-NULL
prob1.est.geom2<-prob2.est.geom2<-prob3.est.geom2<-NULL
prob1.est.geom3<-prob2.est.geom3<-prob3.est.geom3<-NULL
prob1.est.CE<-prob2.est.CE<-prob3.est.CE<-NULL
prob1.est.CE2<-prob2.est.CE2<-prob3.est.CE2<-NULL

thetamat1<-thetamat2<-thetamat3<-matrix(NA,ncol=7,nrow=nrep)

Set1inRegion1<-Set2inRegion1<-Set3inRegion1<-numeric(nrep)

time.start<-Sys.time()
# set.seed(101)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  set.seed(k)
  
  ######## GEOMETRIC ################
  x.norm = qexp(pnorm(rmvnorm(2500,sigma=matrix(c(1,rho,rho,
                                                  rho,1,rho,
                                                  rho,rho,1),3,3))))
  x.alog = rmvevd(2500, dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1))
  x.alog<--log(1-exp(-1/x.alog))
  x = rbind(x.norm,x.alog)
  
  r<-apply(x,1,sum)
  w<-x/r
  
  qr<-QR.3d(r = r, w=w, tau=0.95)
  
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind,]
  # Threshold value corresponding to each w:
  r0w<-qr$r0w[excind]
  na.ind<-which(is.na(rexc))
  na.ind
  if(length(na.ind)>0){
    rexc<-rexc[-na.ind]
    wexc<-wexc[-na.ind,]
    r0w<-r0w[-na.ind]}
  
  # Check if sets are in region R'>1
  Set1inRegion1[k]<-checkset.3d(x=c(x11,x12),y=c(y11,y12),z=c(z11,z12),w=w,r0w=qr$r0w)
  Set2inRegion1[k]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),w=w,r0w=qr$r0w)
  Set3inRegion1[k]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),w=w,r0w=qr$r0w)
  
  # taus
  m31<-method2(x, delta=0.4, hillQuantile=0.9, pi= 0.05)
  m31
  m32<-method2(x, delta=0.5, hillQuantile=0.9, pi= 0.05)
  m32
  m33<-method2(x, delta=0.6, hillQuantile=0.9, pi= 0.05)
  m33
  
  thetavec1<-rep(NA,7)
  thetavec1[m31$prob>0]<-1
  thetavec1[m31$prob==0]<-10e10
  
  thetavec2<-rep(NA,7)
  thetavec2[m32$prob>0]<-1
  thetavec2[m32$prob==0]<-10e10
  
  thetavec3<-rep(NA,7)
  thetavec3[m33$prob>0]<-1
  thetavec3[m33$prob==0]<-10e10
  
  thetamat1[k,]<-thetavec1
  thetamat2[k,]<-thetavec2
  thetamat3[k,]<-thetavec3
  
  skip_to_next <- FALSE
  tryCatch({
    mygfun1<-function(xyz,par)
    {
      gauge_rvad_all3d(xyz=xyz,par=par,theta1=thetavec1[1],theta2=thetavec1[2], theta3=thetavec1[3],theta12 = thetavec1[4],theta13 = thetavec1[5],theta23 = thetavec1[6],theta123 = thetavec1[7])
    }
    
    ind1<-which(thetavec1==1)
    
    f1<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = mygfun1,
                            init.val = c(2,rep(0.5,len=length(ind1[ind1>3]))),
                            lower.limit=rep(0,1+length(ind1[ind1>3])),
                            upper.limit=c(100,rep(1,length(ind1[ind1>3]))),
                            control=list(reltol=1e-6))
    
    if(!identical(thetavec1,thetavec2))
    {
      mygfun2<-function(xyz,par)
      {
        gauge_rvad_all3d(xyz=xyz,par=par,theta1=thetavec2[1],theta2=thetavec2[2],
                         theta3=thetavec2[3],theta12 = thetavec2[4],theta13 = thetavec2[5],
                         theta23 = thetavec2[6],theta123 = thetavec2[7])
      }
      ind2<-which(thetavec2==1)
      
      f2<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = mygfun2,
                              init.val = c(2,rep(0.5,len=length(ind2[ind2>3]))),
                              lower.limit=rep(0,1+length(ind2[ind2>3])),
                              upper.limit=c(100,rep(1,length(ind2[ind2>3]))),
                              control=list(reltol=1e-6))
      } else{
      ind2<-NULL
      f2<-list(mle=NULL,nllh=Inf,convergence=0)
    }
    
    if(!identical(thetavec1,thetavec3)||!identical(thetavec2,thetavec3))
    {
      mygfun3<-function(xyz,par)
      {
        gauge_rvad_all3d(xyz=xyz,par=par,theta1=thetavec3[1],theta2=thetavec3[2],
                         theta3=thetavec3[3],theta12 = thetavec3[4],theta13 = thetavec3[5],
                         theta23 = thetavec3[6],theta123 = thetavec3[7])
      }
      
      ind3<-which(thetavec3==1)
      
      f3<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = mygfun3,
                              init.val = c(2,rep(0.5,len=length(ind3[ind3>3]))),
                              lower.limit=rep(0,1+length(ind3[ind3>3])),
                              upper.limit=c(100,rep(1,length(ind3[ind3>3]))),
                              control=list(reltol=1e-6))
    }else{
      ind3<-NULL
      f3<-list(mle=NULL,nllh=Inf,convergence=0)
    }
    
    mygfun4 = function(xyz,par){
      # gauge for d=3 simulation study
      # par -> length 5
      if(length(par) != 5){
        stop("parameters not of correct length")
      }
      min(geometricMVE::gauge_gaussian3d(xyz,par=par[1:3]),
          geometricMVE::gauge_rvad_all3d(xyz,par=par[4:5],
                                         theta1=0,theta2=0,theta3=0,theta13=0,theta23=0))
    }
    f4 = fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = mygfun4,
                             init.val = c(2,rep(0.5,len=5)),
                             # lower.limit=rep(0,1+length(ind3[ind3>3])),
                             # upper.limit=c(100,rep(1,length(ind3[ind3>3]))),
                             control=list(maxit=1e5))
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    next
  }
  llvals[k,]<-c(f1$nllh,f2$nllh,f3$nllh,f4$nllh)
  AICvals[k,]<-c(2*f1$nllh+2*(1+length(ind1[ind1>3])),2*f2$nllh+2*(1+length(ind2[ind2>3])),2*f3$nllh+2*(1+length(ind3[ind3>3])),2*(f4$nllh + length(f4$mle)))
  
  if(which.min(AICvals[k,])==1){fittedpar<-f1$mle;fittedgfun<-mygfun1}
  if(which.min(AICvals[k,])==2){fittedpar<-f2$mle;fittedgfun<-mygfun2}
  if(which.min(AICvals[k,])==3){fittedpar<-f3$mle;fittedgfun<-mygfun3}
  if(which.min(AICvals[k,])==4){fittedpar<-f4$mle;fittedgfun<-mygfun4}
  
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
  
  xstar3<-sim.3d(k=ks,w = wexc,r0w = r0w,nsim = n1,par = fittedpar, gfun = fittedgfun)
  
  prob1.est.geom3[k]<-mean(xstar3[,1]>x11&xstar3[,1]<x12&xstar3[,2]>y11&xstar3[,2]<y12&xstar3[,3]>z11&xstar3[,3]<z12) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob2.est.geom3[k]<-mean(xstar3[,1]>x21&xstar3[,1]<x22&xstar3[,2]>y21&xstar3[,2]<y22&xstar3[,3]>z21&xstar3[,3]<z22) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  prob3.est.geom3[k]<-mean(xstar3[,1]>x31&xstar3[,1]<x32&xstar3[,2]>y31&xstar3[,2]<y32&xstar3[,3]>z31&xstar3[,3]<z32) * Rexc.prob.k.3d(k=ks, r0w=r0w, w=wexc, gfun=fittedgfun, par=fittedpar)* length(rexc)/length(r)
  
  ###################################
  
  ######## CE #######################
  
  th<-quantile(x[,1],0.95)
  # 2|1
  fce1<-optim(HT.nll.PL,par=c(0.5,0.5),dat=x[,c(1,2)],u=th,which=1)
  Z1<-(x[x[,1]>th,2]-fce1$par[1]*x[x[,1]>th,1])/x[x[,1]>th,1]^fce1$par[2]
  
  #3|1
  fce2<-optim(HT.nll.PL,par=c(0.5,0.5),dat=x[,c(1,3)],u=th,which=1)
  Z2<-(x[x[,1]>th,3]-fce2$par[1]*x[x[,1]>th,1])/x[x[,1]>th,1]^fce2$par[2]
  
  Z<-cbind(Z1,Z2)
  
  n1<-50000
  xe.star<-rexp(n1)+qexp(0.999)
  star.ind<-sample(1:length(Z1),size=n1,replace=T)
  Zstar<-Z[star.ind,]
  Y1.star<-fce1$par[1]*xe.star+Zstar[,1]*xe.star^fce1$par[2]
  Y2.star<-fce2$par[1]*xe.star+Zstar[,2]*xe.star^fce2$par[2]
  xstarce<-cbind(xe.star,Y1.star,Y2.star)
  
  
  prob1.est.CE[k]<-mean(xstarce[,1]>x11&xstarce[,1]<x12&xstarce[,2]>y11&xstarce[,2]<y12&xstarce[,3]>z11&xstarce[,3]<z12) * 0.001
  prob2.est.CE[k]<-mean(xstarce[,1]>x21&xstarce[,1]<x22&xstarce[,2]>y21&xstarce[,2]<y22&xstarce[,3]>z21&xstarce[,3]<z22) * 0.001
  prob3.est.CE[k]<-mean(xstarce[,1]>x31&xstarce[,1]<x32&xstarce[,2]>y31&xstarce[,2]<y32&xstarce[,3]>z31&xstarce[,3]<z32) * 0.001
  
  # Higher threshold
  
  xe.star<-rexp(n1)+8
  star.ind<-sample(1:length(Z1),size=n1,replace=T)
  Zstar<-Z[star.ind,]
  Y1.star<-fce1$par[1]*xe.star+Zstar[,1]*xe.star^fce1$par[2]
  Y2.star<-fce2$par[1]*xe.star+Zstar[,2]*xe.star^fce2$par[2]
  xstarce<-cbind(xe.star,Y1.star,Y2.star)
  
  prob1.est.CE2[k]<-mean(xstarce[,1]>x11&xstarce[,1]<x12&xstarce[,2]>y11&xstarce[,2]<y12&xstarce[,3]>z11&xstarce[,3]<z12) * exp(-8)
  prob2.est.CE2[k]<-mean(xstarce[,1]>x21&xstarce[,1]<x22&xstarce[,2]>y21&xstarce[,2]<y22&xstarce[,3]>z21&xstarce[,3]<z22) * exp(-8)
  prob3.est.CE2[k]<-mean(xstarce[,1]>x31&xstarce[,1]<x32&xstarce[,2]>y31&xstarce[,2]<y32&xstarce[,3]>z31&xstarce[,3]<z32) * exp(-8)
  
  ###################################
  if(k%%20==0){save.image("SimStudy3d_mix.Rdata")}
}
time.end<-Sys.time()

save.image("SimStudy3d_mix.Rdata")