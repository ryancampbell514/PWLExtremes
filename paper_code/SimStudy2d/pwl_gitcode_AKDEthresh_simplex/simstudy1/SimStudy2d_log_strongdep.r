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

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# library(PWLExtremes)
# source("~/GitHub/PWLExtremes/R/likelihoodandmodelfitting.R")
fn.dir = "../../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
source("../checkset.R")

source("../../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")

#----------------------------------

nrep<-200
nnodes=11
par.locs = seq(0,1,length.out=nnodes)
pen.const = 1
pen.const.vals = NULL

pars<-matrix(NA,ncol=nnodes,nrow=nrep)
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

aAD<-0.4
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

load("../../quantiles/quantiles_akde/log_strongdep_quantiles.Rdata")

time.start<-Sys.time()
#set.seed(321)
for(k in 1:nrep){
  message(paste(k,"/",nrep))
  
  set.seed(k)
  x<-log.strongdep.res[[k]]$x
  r<-x[,1]+x[,2]
  w<-x[,1]/r
  
  #==================================
  
  qr<-log.strongdep.res[[k]]$qr
  
  r0w=qr$r0w
  wpts = qr$wpts
  excind<-r>r0w
  rexc<-r[excind]
  wexc<-w[excind]
  r0w<-r0w[excind]
  
  # plot(rexc*cbind(wexc,1-wexc))
  # if(k>5) stop()
  # next
  
  ##################################

  skip_to_next <- FALSE
  tryCatch({
    model.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=pen.const,fixshape=T)  # this is bad!
    # model.fit1 = PWLExtremes::fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=pen.const,method="BFGS")  # this is good!
    # model.fit$mle
    # model.fit1$mle
    # 
    # # check if the gauge functions from geometricMVE and PWLExtremes compute the same values
    # w.adj.angles = which.adj.angles.2d(wexc, par.locs)
    # psi=model.fit1$mle
    # gvals = pwlin.g.vals.2d(w.adj.angles,psi,par.locs)
    # gvals.1 = PWLExtremes::pwlin.g.vals.2d(w.adj.angles,psi,par.locs)
    # # yes
    # 
    # # check if the log-likelihoods are the same
    # nll.pwlin.2d(psi,rexc,r0w,w.adj.angles,par.locs,W.fit=F,joint.fit=F,
    #                        pen.const=1,fixshape=T)
    # PWLExtremes::nll.pwlin.2d(psi,rexc,r0w,w.adj.angles,par.locs,norm, marg = "pos", 
    #                           pen.norm = "2", 
    #                           pen.adj=F, fW.fit=F, joint.fit=F, pos.par = TRUE, pen.const = 1, 
    #                           fixed.shape = TRUE)
    # # yes they are the same
    
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    message("optim error, returning NAs.")
    next
  }
  
  pars[k,]<-model.fit$mle
  nll.vals[k] = model.fit$nll
  pen.const.vals = c(pen.const.vals,model.fit$pen.const)
  
  ##################################

  xstar<-sim.geometric(model.fit,nsim=n1,k=1)
  pwl.prob1.est.geom[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)
  pwl.prob2.est.geom[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * length(rexc)/length(r)
  pwl.prob3.est.geom[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * length(rexc)/length(r)
  
  
  ## Extrapolation from higher levels
  
  # Find appropriate values of k (for R'>k)
  # k.val=3
  wpts=seq(0,1,length.out=100)
  r.tau.wpts = eval.thresh(qr,wpts)
  cm<-max(max((wpts)*r.tau.wpts,na.rm=T),max((1-wpts)*r.tau.wpts,na.rm=T),na.rm=T)
  km<-10/cm #  max value that works everywhere given the regions we're interested in
  pr.k.val = mean(pgamma(km*r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = model.fit$mle), lower.tail = F)/
                pgamma(r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = model.fit$mle), lower.tail = F))
  xstar<-sim.geometric(model.fit,nsim=n1,k=km)
  pwl.prob1.est.geom2[k]<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * pr.k.val* length(rexc)/length(r)
  pwl.prob2.est.geom2[k]<-mean(xstar[,1]>x21&xstar[,1]<x22&xstar[,2]>y21&xstar[,2]<y22) * pr.k.val* length(rexc)/length(r)
  pwl.prob3.est.geom2[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * pr.k.val* length(rexc)/length(r)
  
  # Higher level just for set 3
  # k.val = 4
  cs<-NULL
  kgrid<-seq(km,10,len=100)#seq(km,4.2,len=10)
  for(j in 1:100)
  {
    cs[j]<-checkset.2d(x=c(x31,x32),y=c(y31,y32),r0w = qr$r0w,w=w,k = kgrid[j])
  }
  ks<-max(kgrid[cs],na.rm=T)
  pr.k.val = mean(pgamma(ks*r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = model.fit$mle), lower.tail = F)/
                    pgamma(r0w, shape = model.fit$shape, rate = sapply(wexc, gfun, par = model.fit$mle), lower.tail = F))
  xstar<-sim.geometric(model.fit,nsim=n1,k=ks)
  pwl.prob3.est.geom3[k]<-mean(xstar[,1]>x31&xstar[,1]<x32&xstar[,2]>y31&xstar[,2]<y32) * pr.k.val* length(rexc)/length(r)
  
  
  #########
  if(k%%20==0){save.image("pwlinexp_SimStudy2d_log_strongdep.Rdata")}
}

time.end<-Sys.time()

rm(log.strongdep.res)
save.image("pwlinexp_SimStudy2d_log_strongdep.Rdata")

stop()
