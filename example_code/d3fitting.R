# d=3 example script for model fitting in the piecewise-linear geometric setting

rm(list=ls())

library(geometry)
library(geometricMVE)
library(evd)
library(mvtnorm)
library(rgl)
library(lattice)
library(PWLExtremes)

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T)

fn.dir = "~/GitHub/PWLExtremes/R"  #PATH TO PWLEXTREMES/R
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# setwd("~/GitHub/PWLExtremes/example_code")

set.seed(4444)
n = 5000  # generate n datapoints
tau=0.95  # the quantile at which we estimate the radial threshold

# asymmetric logistic data, pairwise dependence
asy <- list(.0, .0, .0, c(.5,0.5), c(.5,0.5), c(.5,.5), c(.0,.0,.0))
x<-rmvevd(n, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
x<--log(1-exp(-1/x))

# # asymmetric logistic data, pairwise dependence, dependence on (1), (1,2), (2,3)
# asy <- list(0.5, 0, 0, c(0.5,0.5), c(0,0), c(0.5,1), c(.0,.0,.0))
# x<-rmvevd(n, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
# x<--log(1-exp(-1/x))
#
# # Gaussian and asymmetric logistic mixture model
# asy <- list(0, 0, 0, c(0.5,0.5), c(0,0), c(0,0), c(0.5,0.5,1))
# rho=0.4
# phi=0.6
# x.norm = qexp(pnorm(rmvnorm(floor(n/2),sigma=matrix(c(1,rho,rho,
#                                                 rho,1,rho,
#                                                 rho,rho,1),3,3))))
# x.alog = rmvevd(ceiling(n/2), dep = phi, asy = asy, model = "alog", d = 3,mar=c(1,1,1))
# x.alog<--log(1-exp(-1/x.alog))
# x = rbind(x.norm,x.alog)

# obtain radii and angles
r<-apply(x,1,sum)
w<-x/r

# estimate the threshold
qr = radial.quants.L1.KDE(r,w,tau=tau,bww=0.05,bwr=0.05)

# keep the exceedances
excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind,]
r0w<-qr$r0w[excind]
na.ind<-which(is.na(rexc))
if(length(na.ind)>0){
  rexc<-rexc[-na.ind]
  wexc<-wexc[-na.ind,]
  r0w<-r0w[-na.ind]}

# plot the KDE threshold
w.mesh = qr$wpts
n.mesh = sqrt(nrow(w.mesh))
r.tau.mat = matrix(qr$r.tau.wpts, nrow=n.mesh, ncol=n.mesh)
plot3d(x)
surface3d(x=w.mesh[,1]*r.tau.mat,
          y=w.mesh[,2]*r.tau.mat,
          z=w.mesh[,3]*r.tau.mat,
          col="red",alpha=0.5)
axes3d(edges="bbox")
view3d(userMatrix = usermat,zoom=0.8)

# define the reference angles
par.locs = seq(0,1,by=1/6)
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])
# par.locs = get_ref_angles(wexc)

# plot the Delaunay triangulation of the S_2 simplex
par(pty="s",mfrow=c(1,1))#,mar=c(4,4,1,4))
plot(par.locs[,-3],type="p",pch=20,cex=2,xlab=expression(w[1]),ylab=expression(w[2]))
segments(0,0,1,0,lwd=2)
segments(0,0,0,1,lwd=2)
segments(1,0,0,1,lwd=2)
del.tri = geometry::delaunayn(p=par.locs[,-3], output.options=TRUE)
for(i in 1:nrow(del.tri$tri)){
  vertices = del.tri$tri[i,]
  tris = combn(vertices,2)
  apply(tris,2,function(which.pts){
    pts = par.locs[which.pts,-3]
    segments(pts[1,1],pts[1,2],pts[2,1],pts[2,2])
  })
}

# fit the models
model.fit.R.unbounded  = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=0.1,method="BFGS",bound=FALSE)
model.fit.R.bounded    = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=0.1,method="BFGS",bound=TRUE)
model.fit.W            = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=2,fW.fit=T,method="BFGS")
model.fit.RW.unbounded = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=0.1,method="BFGS",bound=FALSE,fW.fit=T,joint.fit=T)
model.fit.RW.bounded   = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,pen.const=0.1,method="BFGS",bound=TRUE,fW.fit=T,joint.fit=T)

# plot the unit level sets
g.vals = gfun.pwl(x=w.mesh,par=model.fit.R.unbounded$mle,ref.angles=par.locs)
g.vals.mat =  matrix(g.vals,n.mesh,n.mesh)
open3d()
plot3d(x/log(n))
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab="x1",ylab="x2",zlab="x3")
view3d(userMatrix = usermat,zoom=0.8)

g.vals = gfun.pwl(x=w.mesh,par=model.fit.R.bounded$mle,ref.angles=par.locs)
g.vals.mat =  matrix(g.vals,n.mesh,n.mesh)
open3d()
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab="x1",ylab="x2",zlab="x3")
view3d(userMatrix = usermat,zoom=0.8)

g.vals = gfun.pwl(x=w.mesh,par=model.fit.RW.unbounded$mle,ref.angles=par.locs)
g.vals.mat =  matrix(g.vals,n.mesh,n.mesh)
open3d()
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab="x1",ylab="x2",zlab="x3")
view3d(userMatrix = usermat,zoom=0.8)

g.vals = gfun.pwl(x=w.mesh,par=model.fit.RW.bounded$mle,ref.angles=par.locs)
g.vals.mat =  matrix(g.vals,n.mesh,n.mesh)
open3d()
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab="x1",ylab="x2",zlab="x3")
view3d(userMatrix = usermat,zoom=0.8)

# plot the fitted angular densities on the S_2 simplex
g.vals.fW = gfun.pwl(x=w.mesh,par=model.fit.W$fW.mle,ref.angles=par.locs)
G.vol.W = G.vol(model.fit.W$fW.mle,par.locs)
d = data.frame(x=w.mesh[,1], y=w.mesh[,2], z=(1/(3*G.vol.W))*(g.vals.fW^(-3)))
par(pty="s",mfrow=c(1,1))
lattice::levelplot(z~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))

g.vals.fW = gfun.pwl(x=w.mesh,par=model.fit.RW.unbounded$mle,ref.angles=par.locs)
G.vol.W = G.vol(model.fit.RW.unbounded$mle,par.locs)
d = data.frame(x=w.mesh[,1], y=w.mesh[,2], z=(1/(3*G.vol.W))*(g.vals.fW^(-3)))
par(pty="s",mfrow=c(1,1))
lattice::levelplot(z~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))

g.vals.fW = gfun.pwl(x=w.mesh,par=model.fit.RW.bounded$mle,ref.angles=par.locs)
G.vol.W = G.vol(model.fit.RW.bounded$mle,par.locs)
d = data.frame(x=w.mesh[,1], y=w.mesh[,2], z=(1/(3*G.vol.W))*(g.vals.fW^(-3)))
par(pty="s",mfrow=c(1,1))
lattice::levelplot(z~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))

# generate some data from 1 of the fitted models
n1 = 50000
gfun = function(w,par,locs=par.locs){
  gfun.pwl(x=w,par=par,ref.angles=locs)
}
xstar = sim.cond(w=wexc, r0w=r0w, nsim=n1,k=1,gfun=gfun,par=model.fit.R.unbounded$mle)
xstar = sim.cond(w=wexc, r0w=r0w, nsim=n1,k=1,gfun=gfun,par=model.fit.R.bounded$mle)
xstar = sim.joint(nsim=n1,k.vals=1,par=model.fit.R.unbounded$mle,fW.par=model.fit.W$fW.mle,par.locs=par.locs,r=r,w=w)[[1]]
xstar = sim.joint(nsim=n1,k.vals=1,par=model.fit.R.bounded$mle,fW.par=model.fit.W$fW.mle,par.locs=par.locs,r=r,w=w)[[1]]
xstar = sim.joint(nsim=n1,k.vals=1,par=model.fit.RW.unbounded$mle,fW.par=model.fit.RW.unbounded$mle,par.locs=par.locs,r=r,w=w)[[1]]
xstar = sim.joint(nsim=n1,k.vals=1,par=model.fit.RW.unbounded$mle,fW.par=model.fit.RW.unbounded$mle,par.locs=par.locs,r=r,w=w)[[1]]

# define a region in which to estimate probabilities
x11<-8
x12<-10
y11<-8
y12<-10
z11<-0.01
z12<-3
prob.est<-mean(xstar[,1]>x11&xstar[,1]<x12&xstar[,2]>y11&xstar[,2]<y12&xstar[,3]>z11&xstar[,3]<z12) * length(rexc)/length(r)

par(mfrow=c(1,1),pty="s")
plot3d(x,xlim=c(0,16),ylim=c(0,16))
points3d(xstar,col="blue")
