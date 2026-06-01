rm(list=ls())

library(evd)
library(qgam)
library(rmutil)
library(geometry)
library(scales)
library(geometricMVE)
library(mvtnorm)
library(rgl)

# # load-in the functions
# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
library(PWLExtremes)

################################################################################

# read-in d=2 data with logistic dependence
nnodes=11
par.locs = seq(0,1,length.out=nnodes)
pen.const = 0.1

set.seed(909)  # best is 1111
x<-rbvevd(5000,dep=0.4,mar1=c(0,1,0))
x<-qexp(evd::pgumbel(x))
r<-x[,1]+x[,2]
w<-x[,1]/r

tau=0.95
qr = PWLExtremes::radial.quants.L1.KDE.2d(r,w,tau=tau,bww=0.05,bwr=0.05,ker="Gaussian")

r0w=qr$r0w
wpts = qr$wpts
wpts[101] =0.5

excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind]
r0w<-r0w[excind]

# Fit the PWL models
init = KDE.quant.eval.2d(wpts=par.locs,r=r,w=w,tau=tau,ker="Gaussian")
init = init/max(par.locs * init, na.rm = TRUE)

## unbounded model, no pen
model.fit.nopen13 = PWLExtremes::fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
                                              pen.const=0,method="BFGS",init.val=init)

## unbounded model
model.fit.pen24 = PWLExtremes::fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
                                            pen.const=1,method="BFGS",init.val=init,bound.fit=TRUE)
# ## bounded model
# model.fit.pen24 = bound.fit.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
#                                pen.const=pen.const,method="BFGS",init.val=init)
# ## angular model
# fW.fit34 = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
#                         pen.const=pen.const,method="BFGS",fW.fit=T,init.val=init)
# ## unbounded joint model
# model.fit.pen5 = fit.pwlin.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
#                               pen.const=pen.const,method="BFGS",
#                               fW.fit=T,joint.fit=T,init.val=init)
# ## bounded joint model
# model.fit.pen6 = bound.fit.2d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,
#                               pen.const=pen.const,method="BFGS",
#                               fW.fit=T,joint.fit=T,init.val=init)

# plot gauges
# par(mfrow=c(1,2),pty="s")
pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs/d2log_nopen.pdf",width=4,height=4)
par(mfrow=c(1,1),pty="s")
plot(x/log(5000),col="grey",pch=20,
     xlab=expression(x[1]/log(n)),
     ylab=expression(x[2]/log(n)),
     xlim=c(0,1.2),ylim=c(0,1.2))
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.nopen13$mle,
                                 ref.angles=par.locs),lwd=2)
segments(0,1,1,1,lty=3)
segments(1,0,1,1,lty=3)
dev.off()

pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs/d2log_pen_bound.pdf",width=4,height=4)
par(mfrow=c(1,1),pty="s")
plot(x/log(5000),col="grey",pch=20,
     xlab=expression(x[1]/log(n)),
     ylab=expression(x[2]/log(n)),
     xlim=c(0,1.2),ylim=c(0,1.2))
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.pen24$mle,
                                 ref.angles=par.locs),lwd=2)
segments(0,1,1,1,lty=3)
segments(1,0,1,1,lty=3)
dev.off()

par(mfrow=c(1,4),pty="s")
plot(x/log(5000),col="grey")
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.pen13$mle,
                                 ref.angles=par.locs))
plot(x/log(5000),col="grey")
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.pen24$mle,
                                 ref.angles=par.locs))
plot(x/log(5000),col="grey")
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.pen5$mle,
                                 ref.angles=par.locs))
plot(x/log(5000),col="grey")
lines(cbind(wpts,1-wpts)/gfun.2d(x=cbind(wpts,1-wpts),
                                 par=model.fit.pen6$mle,
                                 ref.angles=par.locs))

##############################################################################

# d=3 Delauney triangulation on a d=3 dataset
rm(list = ls())

library(evd)
library(geometry)
library(pracma)
library(rgl)
library(lattice)
library(geometricMVE)

setwd("~/Dropbox/phd_research/pw_lin_gauge/pen_score_plots/")

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
library(PWLExtremes)

# the perspective  at which we want to view rgl figures
usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

set.seed(11)
# set.seed(311)
asy <- list(.0, .0, .0, c(.5,0.5), c(.5,0.5), c(.5,.5), c(.0,.0,.0))
ds<-rmvevd(5000, dep = c(.4), asy = asy, model = "alog", d = 3,mar=c(1,1,1))
ds<--log(1-exp(-1/ds))
r<-apply(ds,1,sum)
w<-ds/r

tau=0.95
qr=radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=T,ker="Gaussian")
r0w = qr$r0w
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind,]
r0wexc<-r0w[excind]

par.locs = seq(0,1,by=1/6)   # newly-suggested 
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])
nnodes = nrow(par.locs)

init = KDE.quant.eval(wpts=par.locs,r=r,w=w,tau=tau,ker="Gaussian")
init = init/max(par.locs[,1] * init, na.rm = TRUE)

model.fit.24 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,
                         locs=par.locs,pen.const=0,method="BFGS",
                         init.val=init)
model.fit.pen24 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,
                            locs=par.locs,pen.const=1,  #0.003333333
                            method="BFGS",
                            init.val=init,bound.fit=T)
g.est.vals.mod1 = gfun.pwl(x=qr$wpts,par=model.fit.24$mle,ref.angles=par.locs)
g.est.vals.mod2 = gfun.pwl(x=qr$wpts,par=model.fit.pen24$mle,ref.angles=par.locs)

open3d()
plot3d(ds/log(nrow(ds)),xlab=" ",ylab=" ",zlab=" ")
surface3d(x=qr$wpts[,1]/matrix(g.est.vals.mod1,nrow=30,ncol=30),
          y=qr$wpts[,2]/matrix(g.est.vals.mod1,nrow=30,ncol=30),
          z=qr$wpts[,3]/matrix(g.est.vals.mod1,nrow=30,ncol=30),
          col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat,zoom=0.8)
points3d(1.5*diag(3),cex=0)  # to make x,y,zlims c(0,1.5)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/d3alog1_nopen.png",width=500,height=500)
snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs/d3alog1_nopen.png",width=500,height=500)

open3d()
plot3d(ds/log(nrow(ds)),xlab=" ",ylab=" ",zlab=" ")
surface3d(x=qr$wpts[,1]/matrix(g.est.vals.mod2,nrow=30,ncol=30),
          y=qr$wpts[,2]/matrix(g.est.vals.mod2,nrow=30,ncol=30),
          z=qr$wpts[,3]/matrix(g.est.vals.mod2,nrow=30,ncol=30),
          col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat,zoom=0.8)
points3d(1.5*diag(3),cex=0)  # to make x,y,zlims c(0,1.5)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/d3alog1_pen_bound.png",width=500,height=500)
snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs/d3alog1_pen_bound.png",width=500,height=500)