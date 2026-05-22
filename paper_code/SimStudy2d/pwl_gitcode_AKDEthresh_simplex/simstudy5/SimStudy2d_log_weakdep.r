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
source("../checkset.R")

# load-in some needed scripts that aren't on geometricMVE
source("../../../../geometricMVEdev/thresh_estimation_AKDE.R")
source("../../../../geometricMVEdev/sim_joint.R")

source("../../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")
source("../../../sim_joint.R")

#----------------------------------

nrep<-200
nnodes=11
par.locs = seq(0,1,length.out=nnodes)
pen.const = 1
pen.const.W = 20
pen.const.vals = NULL
pen.const.W.vals = NULL

pars<-matrix(NA,ncol=nnodes,nrow=nrep)
fW.pars<-matrix(NA,ncol=nnodes,nrow=nrep)
nll.vals = rep(NA,nrep)

logistic.cdf <- function(x.vec,dep.par){
  # input - vector in Exp(1) margins
  exp(-(sum((-log(1-exp(-x.vec)))^(1/dep.par)))^dep.par)
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

pwl.prob1.est.geom<-pwl.prob1.est.geom2<-rep(NA,nrep)
pwl.prob2.est.geom<-pwl.prob2.est.geom2<-rep(NA,nrep)
pwl.prob3.est.geom<-pwl.prob3.est.geom2<-pwl.prob3.est.geom3<-rep(NA,nrep)

aAD<-0.8
n<-5000
n1<-50000
ww = seq(0,1,by=0.001)

prob1.true = true.prob.in.B.d2(x11,x12,y11,y12,
                               distn.fn=logistic.cdf,dep.par=aAD)
prob2.true = true.prob.in.B.d2(x21,x22,y21,y22,
                               distn.fn=logistic.cdf,dep.par=aAD)
prob3.true = true.prob.in.B.d2(x31,x32,y31,y32,
                               distn.fn=logistic.cdf,dep.par=aAD)

gfun = function(w,par,locs=par.locs){
  adj.angles = which.adj.angles.2d(angles=w, locs)
  pwlin.g.vals.2d(adj.angles,par,locs)
}

time.start<-Sys.time()
# set.seed(321)
load("../../quantiles/quantiles_kde/log_weakdep_quantiles.Rdata")

time.start<-Sys.time()
#set.seed(321)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  set.seed(k)
  x<-log.weakdep.res[[k]]$x
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  
  #==================================
  
  qr<-log.weakdep.res[[k]]$qr
  r0w = qr$r0w
  wpts = qr$wpts
  
  qr$bww = as.numeric(pdfCluster::kepdf(x = w, kernel = "gaussian", bwtype="adaptive",
                                        alpha=0.5)@par$hx)  # need this here.
  
  excind<-r>r0w
  rexc<-r[excind]
  wexc<-w[excind]
  r0w<-r0w[excind]
  
  # wexc.adj.angles = which.adj.angles.2d(angles=wexc, par.locs)
  
  ##################################
  
  skip_to_next <- FALSE
  tryCatch({
    model.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=pen.const,bound.fit=F,fixshape=T,joint.fit=T,W.fit=T)
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    message("optim error, returning NAs.")
    next
  }
  
  pars[k,]<-model.fit$mle
  fW.pars[k,]=model.fit$mle.W
  nll.vals[k] = model.fit$nll
  pen.const.vals = c(pen.const.vals,model.fit$pen.const)

  mle = model.fit$mle
  mle.fW = model.fit$mle.W
  
  ##################################
  
  # compute k values
  wpts=seq(0,1,length.out=100)
  r.tau.wpts = eval.thresh(qr,wpts)
  cm = max(max((wpts)*r.tau.wpts,na.rm=T),max((1-wpts)*r.tau.wpts,na.rm=T),na.rm=T)  
  cs<-NULL
  kgrid<-seq(10/cm,10,len=100)#seq(km,4.2,len=10)
  for(j in 1:100)
  {
    cs[j]<-checkset.2d(x=c(x31,x32),y=c(y31,y32),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  k.vals = c(1,10/cm,max(kgrid[cs],na.rm=T))
  
  xstar.sims = sim.joint.2d(nsim=n1,k.vals=k.vals,shape=model.fit$shape,
                         par.locs=par.locs,par.locs.W=par.locs,
                         par=mle, fW.par=mle.fW,
                         fitted.thresh=qr)  
  ##################################
  
  
  xstar<-xstar.sims[[1]]
  pwl.prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)
  pwl.prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * length(rexc)/length(r)
  pwl.prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * length(rexc)/length(r)
  
  pr.k.val = mean(pgamma(k.vals[2]*r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = mle), lower.tail = F)/
                    pgamma(r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = mle), lower.tail = F))
  xstar<-xstar.sims[[2]]
  pwl.prob1.est.geom2[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * pr.k.val* length(rexc)/length(r)
  pwl.prob2.est.geom2[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * pr.k.val* length(rexc)/length(r)
  pwl.prob3.est.geom2[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * pr.k.val* length(rexc)/length(r)
  
  pr.k.val = mean(pgamma(k.vals[3]*r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = mle), lower.tail = F)/
                    pgamma(r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = mle), lower.tail = F))
  xstar<-xstar.sims[[3]]
  pwl.prob3.est.geom3[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * pr.k.val* length(rexc)/length(r)
  
  #########
  if(k%%20==0){save.image("pwlinexp_SimStudy2d_log_weakdep.Rdata")}
}

time.end<-Sys.time()

rm(log.weakdep.res)
save.image("pwlinexp_SimStudy2d_log_weakdep.Rdata")
