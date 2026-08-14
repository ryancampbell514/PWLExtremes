# d=2 example script for model fitting in the piecewise-linear geometric setting

rm(list=ls())

library(geometry)
library(geometricMVE)
library(evd)
library(mvtnorm)

source(file.path("PATH TO PWLExtremes","extra-functions.R"))

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
qr = fit.thresh(r,w,tau=tau,bww=0.05)

# keep the exceedances
r0w=qr$r0w
wpts = qr$wpts
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind]
r0w<-r0w[excind]

# plot the KDE threshold
par(mfrow=c(1,1),pty="s")
plotfittedthresh(qr)

# Fit the models
par.locs = seq(0,1,length.out=11)
par.locs = seq(0,1,length.out=11)
model.fit.R.unbounded           = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",bound.fit=F,fixshape = T)
model.fit.R.unbounded2          = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",bound.fit=F,fixshape = F)
model.fit.R.bounded             = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",bound.fit=T,fixshape = T)
model.fit.R.bounded2            = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",bound.fit=T,fixshape = F)
model.fit.RW.unbounded          = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",W.fit=T,joint.fit=T,fixshape = T)
model.fit.RW.unbounded2         = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",W.fit=T,joint.fit=T,fixshape = F)
model.fit.RW.bounded            = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",W.fit=T,joint.fit=T,bound.fit=T,fixshape = T)
model.fit.RW.bounded2           = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=1,method="BFGS",W.fit=T,joint.fit=T,bound.fit=T,fixshape = F)
model.fit.W                     = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,thresh.fit=qr,locs=par.locs,pen.const=20,W.fit=T,method="BFGS")

# plot the unit level sets
par(mfrow=c(3,3),pty="s")
plotfittedgauge(model.fit.R.unbounded)
plotfittedgauge(model.fit.R.unbounded2)
plotfittedgauge(model.fit.R.bounded)
plotfittedgauge(model.fit.R.bounded2)
plotfittedgauge(model.fit.RW.unbounded)
plotfittedgauge(model.fit.RW.unbounded2)
plotfittedgauge(model.fit.RW.bounded)
plotfittedgauge(model.fit.RW.bounded2)

# plot the angular models
wpts = seq(0,1,length.out=100)
par(mfrow=c(1,3),pty="s")
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.W$mle.W,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.W$mle.W, par.locs = par.locs)))
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.unbounded$mle,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.RW.unbounded$mle, par.locs = par.locs)))
hist(wexc,main=NULL,freq=F,ylim=c(0,5))
lines(wpts, gfun.2d(cbind(wpts,1-wpts),par=model.fit.RW.bounded$mle,ref.angles=par.locs)^(-2) / (2*G.vol.2d(gauge.pars=model.fit.RW.bounded$mle, par.locs = par.locs)))

# generate some data from 1 of the fitted models
n1 = 50000
xstar = sim.geometric(fit=model.fit.R.unbounded,nsim=n1,k=1,)

# define a region in which to estimate probabilities
x11<-10
x12<-12
y11<-10
y12<-12
prob.est<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12) * length(rexc)/length(r)

par(mfrow=c(1,1),pty="s")
plot(x,xlim=c(0,16),ylim=c(0,16))
points(xstar,col="blue")
