rm(list=ls())

library(evd)
library(qgam)
library(rmutil)
library(geometry)
library(scales)
# library(geometricMVE)
library(mvtnorm)
library(this.path)

par(pty="s")
setwd(this.path::here())

#----------------------------------

# load-in my version of geometricMVE
fn.dir = "../../../../geometricMVE_RCcode/"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# load-in some needed scripts that aren't on geometricMVE
source("../../../../geometricMVEdev/checkset.R")
source("../../../../geometricMVEdev/thresh_estimation_AKDE.R")
source("../../../../geometricMVEdev/sim_joint.R")

#----------------------------------

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

llvals<-AICvals<-matrix(NA,ncol=3,nrow=nrep)
pwl.prob1.est.geom<-pwl.prob2.est.geom<-pwl.prob3.est.geom<-rep(NA,nrep)
pwl.prob1.est.geom2<-pwl.prob2.est.geom2<-pwl.prob3.est.geom2<-rep(NA,nrep)
pwl.prob1.est.geom3<-pwl.prob2.est.geom3<-pwl.prob3.est.geom3<-rep(NA,nrep)

par.locs = seq(0,1,by=1/6)   # newly-suggested 
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])

gfun = function(w,par,locs=par.locs){
  adj.angles = which.adj.angles(angles=w, locs)
  pwlin.g.vals(adj.angles,par,locs)
}

pars = matrix(data=NA,nrow=200,ncol=nrow(par.locs))

# pen.const=0.01333333 #0.005
pen.const.vals = NULL

load("~/lancaster/pw_lin_gauge/SimStudy_3d/quantiles/quantiles_kde/mix_quantiles.Rdata")
time.start<-Sys.time()
# set.seed(101)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  set.seed(k)

  x = mix.res[[k]]$x
  r<-apply(x,1,sum)
  w<-x/r
  
  tau=0.95  # 0.9
  qr = mix.res[[k]]$qr
  
  qr$bww = pdfCluster::kepdf(x = w[,-3], kernel = "gaussian", bwtype="adaptive",
                             alpha=0.5)@par$hx  # need this here.
  
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind,]
  # Threshold value corresponding to each w:
  r0w<-qr$r0w[excind]
  na.ind<-which(is.na(rexc))
  # na.ind
  if(length(na.ind)>0){
    rexc<-rexc[-na.ind]
    wexc<-wexc[-na.ind,]
    r0w<-r0w[-na.ind]}

  ###################################################################
  
  skip_to_next <- FALSE
  tryCatch({
    # model.fit = fit.pwlin(r=rexc,r0w=r0w,w=wexc,
    #                          locs=par.locs,pen.const=1,method="BFGS")
    model.fit = fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,locs=par.locs,pen.const=1)
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    message("optim error, returning NAs.")
    pars[k,]<-rep(NA,length(pars[k,]))
    pwl.prob1.est.geom[k]=NA
    pwl.prob2.est.geom[k]=NA
    pwl.prob3.est.geom[k]=NA
    pwl.prob1.est.geom2[k]=NA
    pwl.prob2.est.geom2[k]=NA
    pwl.prob3.est.geom2[k]=NA
    pwl.prob1.est.geom3[k]=NA
    pwl.prob2.est.geom3[k]=NA
    pwl.prob3.est.geom3[k]=NA
    next
  }
  mle=model.fit$mle
  pars[k,] = mle
  
  ##################################################################
  
  # xstar<-sim.cond(w = wexc,r0w = r0w,nsim = n1,par = mle, gfun = gfun)
  xstar = sim.geometric(fit=model.fit,nsim=n1,k=1)
  pwl.prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12&xstar[,3]>z11&xstar[,3]<z12) * length(rexc)/length(r)
  pwl.prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22&xstar[,3]>z21&xstar[,3]<z22) * length(rexc)/length(r)
  pwl.prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32&xstar[,3]>z31&xstar[,3]<z32) * length(rexc)/length(r)
  
  # Simulate above higher thresholds (for sets 2 and 3)
  # cm<-max(max(c(qr$wpts[,1]*qr$r.tau.wpts),na.rm=T),max(qr$wpts[,2]*qr$r.tau.wpts,na.rm=T),max((1-qr$wpts[,1]-qr$wpts[,2])*qr$r.tau.wpts, na.rm=T))
  wpts=seq(0,1,length.out=50)
  wpts = expand.grid(wpts,wpts)
  wpts = cbind(wpts,1-apply(wpts,1,sum))
  is.valid = wpts[,3]>=0
  wpts=wpts[is.valid,]
  r.tau.wpts = eval.thresh(qr,wpts)
  cm = max(max((wpts[,1])*r.tau.wpts,na.rm=T),
           max((wpts[,2])*r.tau.wpts,na.rm=T),
           max((wpts[,3])*r.tau.wpts,na.rm=T),na.rm=T)  
  km<-8/cm #  max value that works everywhere given the regions we're interested in
  
  # Higher level for set 2
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  ks = ifelse(ks<1,1,ks)
  
  # xstar2<-sim.cond(k=ks, w = wexc,r0w = r0w,nsim = n1,par = mle, gfun = gfun)
  xstar2 = sim.geometric(fit=model.fit,nsim=n1,k=ks)
  
  pr.ks = mean(PWLExtremes::iweights.pwl(k=ks,r0w=r0w,w=wexc,gfun=gfun,par=mle,shape=dim(w)[2]),na.rm=T)
  pwl.prob1.est.geom2[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12&xstar2[,3]>z11&xstar2[,3]<z12) * pr.ks* length(rexc)/length(r)
  pwl.prob2.est.geom2[k]<-mean(xstar2[,1]>x21&xstar2[,1]<x22&xstar2[,2]>y21&xstar2[,2]<y22&xstar2[,3]>z21&xstar2[,3]<z22) * pr.ks* length(rexc)/length(r)
  pwl.prob3.est.geom2[k]<-mean(xstar2[,1]>x31&xstar2[,1]<x32&xstar2[,2]>y31&xstar2[,2]<y32&xstar2[,3]>z31&xstar2[,3]<z32) * pr.ks* length(rexc)/length(r)
  
  
  # Higher level just for set 3
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  ks = ifelse(ks<1,1,ks)
  
  # xstar3<-sim.cond(k=ks, w = wexc,r0w = r0w,nsim = n1,par = mle, gfun = gfun)
  xstar3 = sim.geometric(fit=model.fit,nsim=n1,k=ks)
  
  pr.ks = mean(PWLExtremes::iweights.pwl(k=ks,r0w=r0w,w=wexc,gfun=gfun,par=mle ,shape=dim(w)[2]),na.rm=T)
  pwl.prob1.est.geom3[k]<-mean(xstar3[,1]>x11&xstar3[,1]<x12&xstar3[,2]>y11&xstar3[,2]<y12&xstar3[,3]>z11&xstar3[,3]<z12) * pr.ks* length(rexc)/length(r)
  pwl.prob2.est.geom3[k]<-mean(xstar3[,1]>x21&xstar3[,1]<x22&xstar3[,2]>y21&xstar3[,2]<y22&xstar3[,3]>z21&xstar3[,3]<z22) * pr.ks* length(rexc)/length(r)
  pwl.prob3.est.geom3[k]<-mean(xstar3[,1]>x31&xstar3[,1]<x32&xstar3[,2]>y31&xstar3[,2]<y32&xstar3[,3]>z31&xstar3[,3]<z32) * pr.ks* length(rexc)/length(r)
  
  ###################################
  if(k%%20==0){save.image("SimStudy3d_mix.Rdata")}
}
time.end<-Sys.time()
rm(mix.res)
save.image("SimStudy3d_mix.Rdata")
