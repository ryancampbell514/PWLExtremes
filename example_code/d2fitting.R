# d=2 example script for model fitting in the piecewise-linear geometric setting

rm(list=ls())

library(geometry)
library(geometricMVE)
library(evd)
library(mvtnorm)

library(PWLExtremes)

set.seed(4444)
n = 5000  # generate n datapoints
tau=0.95  # the quantile at which we estimate the radial threshold

# Generate data in exponential margins

# logsitic data
x<-evd::rbvevd(n,dep=0.4,mar1=c(0,1,0))
x<-qexp(evd::pgumbel(x))

# # Gaussian data
# x<-qexp(pnorm(mvtnorm::rmvnorm(5000,sigma=matrix(c(1,0.8,0.8,1),2,2))))
#
# # inverted logistic data
# x<-rbvevd(5000,dep=0.7,mar1=c(1,1,1))
# x<-1/x

# obtain radii and angles
r<-x[,1]+x[,2]
w<-x[,1]/r

# estimate the threshold
qr = radial.quants.L1.KDE.2d(r,w,tau=tau,bww=0.05,bwr=0.05,ker="Gaussian")

# keep the exceedances
r0w=qr$r0w
wpts = qr$wpts
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind]
r0w<-r0w[excind]

# plot the KDE threshold
par(mfrow=c(1,1),pty="s")
plot(x,pch=20,col="grey",xlim=c(0,11),ylim=c(0,11))
lines(cbind(wpts,1-wpts) * qr$r.tau.wpts,lwd=2,col="red")

# Fit the models
par.locs = seq(0,1,length.out=11)
model.fit.R.unbounded           = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",bound.fit=F)
model.fit.R.unbounded2          = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",bound.fit=F,fixed.shape = F)
model.fit.R.unbounded.pensearch = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=NULL,method="BFGS",bound.fit=F)
model.fit.R.bounded             = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",bound.fit=T)
model.fit.R.bounded2            = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",bound.fit=T,fixed.shape = F)
model.fit.RW.unbounded          = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",fW.fit=T,joint.fit=T)
model.fit.RW.unbounded2         = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",fW.fit=T,joint.fit=T,fixed.shape = F)
model.fit.RW.bounded            = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",fW.fit=T,joint.fit=T,bound.fit=T)
model.fit.RW.bounded2           = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=1,method="BFGS",fW.fit=T,joint.fit=T,bound.fit=T,fixed.shape = F)
model.fit.W                     = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=NULL,fW.fit=T,method="BFGS")

# plot the unit level sets
par(mfrow=c(2,2),pty="s")
plot(cbind(wpts,1-wpts)/gfun.2d(cbind(wpts,1-wpts),par=model.fit.R.unbounded$mle,ref.angles=par.locs),type="l")
plot(cbind(wpts,1-wpts)/gfun.2d(cbind(wpts,1-wpts),par=model.fit.R.unbounded.pensearch$mle,ref.angles=par.locs),type="l")
plot(cbind(wpts,1-wpts)/gfun.2d(cbind(wpts,1-wpts),par=model.fit.R.bounded$mle,ref.angles=par.locs),type="l")
plot(cbind(wpts,1-wpts)/gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.unbounded$mle,ref.angles=par.locs),type="l")
plot(cbind(wpts,1-wpts)/gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.bounded$mle,ref.angles=par.locs),type="l")

# plot the angular models
par(mfrow=c(1,3),pty="s")
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.W$fW.mle,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.W$fW.mle, par.locs = par.locs)))
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.unbounded$mle,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.RW.unbounded$mle, par.locs = par.locs)))
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.bounded$mle,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.RW.bounded$mle, par.locs = par.locs)))

# generate some data from 1 of the fitted models
n1 = 50000
gfun = function(w,par,locs=par.locs){
  gfun.2d(x=w,par=par,ref.angles=locs)
}
xstar = sim.2d.cond(w=wexc, r0w=r0w, nsim=n1,k=1,gfun=gfun,par=model.fit.R.unbounded$mle)
# xstar = sim.2d.cond(w=wexc, r0w=r0w, nsim=n1,k=1,gfun=gfun,par=model.fit.R.bounded$mle)
# xstar = sim.2d.joint(nsim=n1,k.vals=1,gfun=gfun,par=model.fit.R.unbounded$mle,fW.par=model.fit.W$fW.mle,par.locs=par.locs,r=r,w=w)[[1]]
# xstar = sim.2d.joint(nsim=n1,k.vals=1,gfun=gfun,par=model.fit.R.bounded$mle,fW.par=model.fit.W$fW.mle,par.locs=par.locs,r=r,w=w)[[1]]
# xstar = sim.2d.joint(nsim=n1,k.vals=1,gfun=gfun,par=model.fit.RW.unbounded$mle,fW.par=model.fit.RW.unbounded$mle,par.locs=par.locs,r=r,w=w)[[1]]
# xstar = sim.2d.joint(nsim=n1,k.vals=1,gfun=gfun,par=model.fit.RW.unbounded$mle,fW.par=model.fit.RW.unbounded$mle,par.locs=par.locs,r=r,w=w)[[1]]

# define a region in which to estimate probabilities
x11<-10
x12<-12
y11<-10
y12<-12
prob.est<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)

par(mfrow=c(1,1),pty="s")
plot(x,xlim=c(0,16),ylim=c(0,16))
points(xstar,col="blue")
