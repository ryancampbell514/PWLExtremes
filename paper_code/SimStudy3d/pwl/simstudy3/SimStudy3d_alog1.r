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

llvals<-matrix(NA,ncol=3,nrow=nrep)
pwl.prob1.est.geom<-pwl.prob2.est.geom<-pwl.prob3.est.geom<-rep(NA,nrep)
pwl.prob1.est.geom2<-pwl.prob2.est.geom2<-pwl.prob3.est.geom2<-rep(NA,nrep)
pwl.prob1.est.geom3<-pwl.prob2.est.geom3<-pwl.prob3.est.geom3<-rep(NA,nrep)

par.locs = seq(0,1,by=1/6)   # newly-suggested 
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])

# par.locs.W = rbind(diag(3),rep(1/3,3),
#                    c(0.5,0.5,0),
#                    c(0.5,0,0.5),
#                    c(0,0.5,0.5))
par.locs.W = par.locs

gfun = function(w,par,locs=par.locs){
  adj.angles = which.adj.angles(angles=w, locs)
  pwlin.g.vals(adj.angles,par,locs)
}

pars = matrix(data=NA,nrow=200,ncol=nrow(par.locs))
fW.pars = matrix(data=NA,nrow=200,ncol=nrow(par.locs))

# pen.const=0.003333333 #0.005
# pen.const.W = 0.2
pen.const.vals = NULL
pen.const.W.vals = NULL


load("~/lancaster/pw_lin_gauge/SimStudy_3d/quantiles/quantiles_kde/alog1_quantiles.Rdata")
load("pwlinexp_SimStudy3d_alog1.Rdata")
start.k = k
time.start<-Sys.time()
# set.seed(101)
for(k in start.k:nrep){
  message(paste(k,"/",nrep))
  set.seed(k)
  
  x<-alog1.res[[k]]$x
  r<-apply(x,1,sum)
  w<-x/r
  
  tau=0.95  # 0.9
  qr = alog1.res[[k]]$qr
  qr$bww = pdfCluster::kepdf(x = w[,-3], kernel = "gaussian", bwtype="adaptive",
                             alpha=0.5)@par$hx  # need this here.
  excind<-r>qr$r0w
  rexc<-r[excind]
  wexc<-w[excind,]
  r0w<-qr$r0w[excind]
  na.ind<-which(is.na(rexc))
  if(length(na.ind)>0){
    rexc<-rexc[-na.ind]
    wexc<-wexc[-na.ind,]
    r0w<-r0w[-na.ind]}

  ###################################################################
  
  skip_to_next <- FALSE
  tryCatch({
    fW.fit = fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,locs=par.locs.W,pen.const=20,W.fit=T)
    model.fit = fit.geometric.pwl(r=rexc,w=wexc,r0w=r0w,locs=par.locs,pen.const=1)
  },error=function(e){
    skip_to_next <<- TRUE 
  })
  if(skip_to_next) {
    next
  }
  
  # print(opt)
  mle=model.fit$mle
  mle.fW = fW.fit$mle.W
  pars[k,] = mle
  fW.pars[k,] = mle.fW
  
  pen.const.vals = c(pen.const.vals,model.fit$pen.const)
  pen.const.W.vals = c(pen.const.W.vals,fW.fit$pen.const)
  
  ###################################################################
  
  # # Check if sets are in regions R'>1, R'>2
  # Set1inRegion1[k]<-checkset.3d(x=c(x11,x12),y=c(y11,y12),z=c(z11,z12),w=w,r0w=qr$r0w)
  # Set2inRegion1[k]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),w=w,r0w=qr$r0w)
  # Set3inRegion1[k]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),w=w,r0w=qr$r0w)

  # get k values, do simulation
  wpts=seq(0,1,length.out=50)
  wpts = expand.grid(wpts,wpts)
  wpts = cbind(wpts,1-apply(wpts,1,sum))
  is.valid = wpts[,3]>=0
  wpts=wpts[is.valid,]
  r.tau.wpts = eval.thresh(qr,wpts)
  cm = max(max((wpts[,1])*r.tau.wpts,na.rm=T),
           max((wpts[,2])*r.tau.wpts,na.rm=T),
           max((wpts[,3])*r.tau.wpts,na.rm=T),na.rm=T)  
  km<-8/cm
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x21,x22),y=c(y21,y22),z=c(z21,z22),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  k.vals = c(1,ks)
  
  cs<-NULL
  kgrid<-seq(km,2.5,len=10)
  for(j in 1:10)
  {
    cs[j]<-checkset.3d(x=c(x31,x32),y=c(y31,y32),z=c(z31,z32),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs])
  
  k.vals = c(k.vals,ks)
  k.vals = ifelse(k.vals<1,1,k.vals)
  
  xstar.sims = sim.joint(nsim=n1,k.vals=k.vals,shape=model.fit$shape,
                            par.locs=par.locs,par.locs.W=par.locs,
                            par=mle, fW.par=mle.fW,
                            fitted.thresh=qr)  
  # xstar.sims = sim.joint(nsim=n1,k.vals=k.vals,par=mle,fW.par=mle.fW,
  #                        par.locs=par.locs,par.locs.W=par.locs.W,r=r,w=w,tau=tau)
  
  ####################################################################
  
  xstar<-xstar.sims[[1]]
  pwl.prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12&xstar[,3]>z11&xstar[,3]<z12) * length(rexc)/length(r)
  pwl.prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22&xstar[,3]>z21&xstar[,3]<z22) * length(rexc)/length(r)
  pwl.prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32&xstar[,3]>z31&xstar[,3]<z32) * length(rexc)/length(r)
  
  # Higher level for set 2
  xstar2<-xstar.sims[[2]]
  pr.ks = mean(PWLExtremes::iweights.pwl(k=k.vals[2],r0w=r0w,w=wexc,gfun=gfun,par=mle,shape=dim(w)[2]),na.rm=T)
  pwl.prob1.est.geom2[k]<-mean(xstar2[,1]>x11&xstar2[,1]<x12&xstar2[,2]>y11&xstar2[,2]<y12&xstar2[,3]>z11&xstar2[,3]<z12) * pr.ks* length(rexc)/length(r)
  pwl.prob2.est.geom2[k]<-mean(xstar2[,1]>x21&xstar2[,1]<x22&xstar2[,2]>y21&xstar2[,2]<y22&xstar2[,3]>z21&xstar2[,3]<z22) * pr.ks* length(rexc)/length(r)
  pwl.prob3.est.geom2[k]<-mean(xstar2[,1]>x31&xstar2[,1]<x32&xstar2[,2]>y31&xstar2[,2]<y32&xstar2[,3]>z31&xstar2[,3]<z32) * pr.ks* length(rexc)/length(r)
  
  # Higher level just for set 3
  xstar3<-xstar.sims[[3]]
  pr.ks = mean(PWLExtremes::iweights.pwl(k=k.vals[3],r0w=r0w,w=wexc,gfun=gfun,par=mle,shape=dim(w)[2]),na.rm=T)
  pwl.prob1.est.geom3[k]<-mean(xstar3[,1]>x11&xstar3[,1]<x12&xstar3[,2]>y11&xstar3[,2]<y12&xstar3[,3]>z11&xstar3[,3]<z12) * pr.ks* length(rexc)/length(r)
  pwl.prob2.est.geom3[k]<-mean(xstar3[,1]>x21&xstar3[,1]<x22&xstar3[,2]>y21&xstar3[,2]<y22&xstar3[,3]>z21&xstar3[,3]<z22) * pr.ks* length(rexc)/length(r)
  pwl.prob3.est.geom3[k]<-mean(xstar3[,1]>x31&xstar3[,1]<x32&xstar3[,2]>y31&xstar3[,2]<y32&xstar3[,3]>z31&xstar3[,3]<z32) * pr.ks* length(rexc)/length(r)
  
  ###################################
  if(k%%20==0){save.image("pwlinexp_SimStudy3d_alog1.Rdata")}
}
time.end<-Sys.time()
rm(alog1.res)
save.image("pwlinexp_SimStudy3d_alog1.Rdata")
