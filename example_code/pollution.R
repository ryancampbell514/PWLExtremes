# d=4 model fitting on air pollution data
rm(list = ls())

library(geometry)
library(rgl)
library(openair)
library(PWLExtremes)

# fn.dir = "~/GitHub/PWLExtremes/R"  #PATH TO PWLEXTREMES/R
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

setwd("~/GitHub/PWLExtremes/example_code")

load("~/GitHub/PWLExtremes/ds_4d_urban.Rdata")
head(ds.exp.4d)

r<-apply(ds.exp.4d,1,sum)
w<-ds.exp.4d/r

tau=0.70
qr=radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=T)
r0w = qr$r0w
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind,]
r0w<-r0w[excind]

# Do a 2-step prodecure:
# First, fit the model on a sparse set of reference angles.
# Then, use these MLEs as starting values for model fitting on a more
# refined set of reference angles.
par.locs = rbind(diag(4),
                 rep(1/4,4),
                 c(1/3,1/3,1/3,0),c(1/3,1/3,0,1/3),c(1/3,0,1/3,1/3),c(0,1/3,1/3,1/3),
                 c(0.5,0.5,0,0),c(0.5,0,0.5,0),c(0.5,0,0,0.5),c(0,0.5,0.5,0),c(0,0.5,0,0.5),c(0,0,0.5,0.5))

# set the initial parameters of model fitting
init = KDE.quant.eval(wpts=par.locs,r=r,w=w,tau=tau)
init = init/max(par.locs[,1] * init, na.rm = TRUE)
init = ifelse(init > 1 / apply(par.locs,1,max), 1/ apply(par.locs,1,max), init)

R.fit.13 = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,pen.const=0.0005,#0.0002343750,
                        method="BFGS")

tri.mat = geometry::delaunayn(p=par.locs[,-4], output.options=TRUE)$tri
new.mesh = t(apply(tri.mat,1,function(row.idx){
  return(apply(par.locs[row.idx,],2,mean))
}))
new.mesh = new.mesh + rnorm(prod(dim(new.mesh)),sd=0.001)
new.mesh = new.mesh / apply(new.mesh,1,sum)
new.mesh = new.mesh[apply(new.mesh,1,function(w) !any(w<0)),]
new.mesh = new.mesh[apply(new.mesh,1,function(vec) sum(vec>0)>2),]
par.locs.2 = data.frame(rbind(par.locs,new.mesh))
par.locs.2 = as.matrix(par.locs.2[!duplicated(par.locs.2), ])

init.2.mod13 = 1/gfun.pwl(x = par.locs.2, par = R.fit.13$mle, ref.angles = par.locs)
init.2.mod13 = ifelse(init.2.mod13 > 1 / apply(par.locs.2,1,max), apply(par.locs.2,1,max), init.2.mod13)

R.fit.13.refit = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs.2,init.val=init.2.mod13,pen.const=0.0001,method="BFGS")

# save the MLEs and their locations
mle = R.fit.13.refit$mle
par.locs = par.locs.2

#############################################################################

# Plot the gauge function projections

proj.g.fn = function(gfun,w,which.w,nm,...){
  # gfun -> gauge that takes in 4-dim vectors
  # w -> 3-min input
  # which.w -> which index to take min over

  w.inp = matrix(NA,nrow=nm,ncol=4)
  w.inp[,-which.w] = matrix(as.numeric(w),ncol=3,nrow=nm,byrow=T)
  w.inp[,which.w] = seq(0,1.3,length.out=nm)

  return(min(gfun(w.inp,...)))
}

nmesh=50
wpts<-expand.grid(seq(0,1,len=nmesh),seq(0,1,len=nmesh))
wpts = cbind(wpts,1-apply(wpts,1,sum))

t1 = Sys.time()
gvals.est1 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=1,nm=nmesh,par=mle,ref.angles=par.locs)
gvals.est1.mat = matrix(gvals.est1,nmesh,nmesh)
t2 = Sys.time()
print(t2-t1)

gvals.est2 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=2,nm=nmesh,par=mle,ref.angles=par.locs)
gvals.est2.mat = matrix(gvals.est2,nmesh,nmesh)

gvals.est3 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=3,nm=nmesh,par=mle,ref.angles=par.locs)
gvals.est3.mat = matrix(gvals.est3,nmesh,nmesh)

gvals.est4 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=4,nm=nmesh,par=mle,ref.angles=par.locs)
gvals.est4.mat = matrix(gvals.est4,nmesh,nmesh)

denom = log(nrow(ds.exp.4d))
ds.scaled = ds.exp.4d / denom

open3d()
plot3d(ds.scaled[,c(2,3,4)])
surface3d(wpts[,1]/gvals.est1.mat,
          wpts[,2]/gvals.est1.mat,
          wpts[,3]/gvals.est1.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)

open3d()
plot3d(ds.scaled[,c(1,3,4)])
surface3d(wpts[,1]/gvals.est2.mat,
          wpts[,2]/gvals.est2.mat,
          wpts[,3]/gvals.est2.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)

open3d()
plot3d(ds.scaled[,c(1,2,4)])
surface3d(wpts[,1]/gvals.est3.mat,
          wpts[,2]/gvals.est3.mat,
          wpts[,3]/gvals.est3.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)

open3d()
plot3d(ds.scaled[,c(1,2,3)])
surface3d(wpts[,1]/gvals.est4.mat,
          wpts[,2]/gvals.est4.mat,
          wpts[,3]/gvals.est4.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)

#############################################################################

# PP, QQ, and return-level plots

mle=R.fit.13.refit$mle
dirnm="mod1"
# mle = res$RW.fit.5$mle
# par.locs = res$par.locs
par.locs = R.fit.13.refit$par.locs
gfun = function(w,par,locs=par.locs){
  gfun.pwl(x=w,par=par,ref.angles=locs)
}

par(mfrow=c(1,2),pty="s")
geometricMVE::ppdiag.4d(r=rexc,w=wexc,r0w=r0w,par=c(4,mle),gfun=gfun)
geometricMVE::qqdiag.4d(r=rexc,w=wexc,r0w=r0w,par=c(4,mle),gfun=gfun)

# return level plots, evaluated on the dataset's angles
gauge.est.vals.w = gfun.pwl(w,par=mle,ref.angles=par.locs)
T.vals = seq(10,1000,by=10)
est.prop.exc = sapply(T.vals,function(T.val){
  return.radii = qgamma(1-(1/T.val), shape=4, rate=gauge.est.vals.w)
  return(mean(r>return.radii))
})
x.vals = log(T.vals)
y.vals = log(1/est.prop.exc)

load("d4pollution_chi_bs_pwl_tau070.Rdata")
bs.estimates = do.call(rbind,lapply(bs.res.pwl, function(lst) lst$est.log.T.bs))
pwl.CI.l = apply(bs.estimates,2,quantile,probs=0.025,na.rm=T)
pwl.CI.u = apply(bs.estimates,2,quantile,probs=0.975,na.rm=T)

pdf(paste0("figs/",dirnm,"/d4pollution_return.pdf"),width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(x=x.vals,y=y.vals,pch=20,
     xlim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     ylim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     xlab="log(T), true",ylab="log(T), estimated",type="p")
segments(-1,-1,10000,10000)
lines(x.vals,pwl.CI.l,lty=2)
lines(x.vals,pwl.CI.u,lty=2)
dev.off()
