# d=3 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
library(geometricMVE)

library(this.path)
setwd(this.path::here())

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
library(PWLExtremes)

# the perspective  at which we want to view rgl figures
usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

###########################

load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")

r<-apply(ds.exp.3d,1,sum)
w<-ds.exp.3d/r

tau=0.95
# source("~/GitHub/PWLExtremes/R/quantileregression.R")
#qr=PWLExtremes::radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=FALSE,bww=0.0222222222222222,ker="Epanechnikov")
qr=PWLExtremes::radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=TRUE,bww=0.075,ker="Gaussian")
# qr=radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=T,bww=NULL)
# source("~/GitHub/PWLExtremes/R/quantileregression.R")
# qr=radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=T,bww=0.1,
#                         ker="Epanechnikov")
# qr = geometricMVE::QR.3d(r=r,w=w,tau=tau)
r0w = qr$r0w
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind,]
r0wexc<-r0w[excind]
if(any(is.na(rexc))){
 cond = is.na(rexc)
 rexc<-rexc[!cond]
 wexc<-wexc[!cond,]
 r0wexc<-r0wexc[!cond]
}
#save.image(file="d3pollution_fits_tau95_bww_penconst_search.RData")

# load("d3pollution_fits_tau95_bww_penconst_search.RData")

par.locs = seq(0,1,by=1/6)   # newly-suggested 
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])
nnodes = nrow(par.locs)

# locs=par.locs
# which.adj.angles.res = which.adj.angles(wexc,par.locs)
# all.regions = lapply(which.adj.angles.res, function(lst) sort(lst$loc.idx))
# idx.locs = unique(all.regions)
# res = lapply(idx.locs, function(loc) {
#   freq = sum(sapply(which.adj.angles.res, function(lst) all(sort(lst$loc.idx) == loc)))
#   # which.w = sapply(which.adj.angles.res, function(lst) all(sort(lst$loc.idx) == loc))
#   which.w = NULL
#   return(list(reg = locs[loc, ], loc = loc, which.w = which.w, freq = freq))
# })
# wexc[res[[31]]$which.w,]
# n.DT.res = sapply(n.DT.region(which.adj.angles(wexc,par.locs),par.locs), function(lst) lst$freq)
# mean(n.DT.res==1)

init = KDE.quant.eval(wpts=par.locs,r=r,w=w,tau=tau,bww=qr$bww,ker="Gaussian")
init = init/max(par.locs[,1] * init, na.rm = TRUE)
# init=NULL

# cbind(par.locs,1/apply(par.locs,1,max),init,init>(1/apply(par.locs,1,max)))

message("Fitting R|W bounded")
# R.fit.24 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,locs=par.locs,pen.const=NULL,method="BFGS",init.val=init,bound.fit=TRUE)
R.fit.24 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,locs=par.locs,pen.const=1,method="BFGS",init.val=init,bound.fit=TRUE)
# save.image(file="d3pollution_fits_tau95_bww_penconst_search.RData")

# message("Fitting W")
# W.fit.34 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,locs=par.locs,pen.const=NULL,method="BFGS",joint.fit=F,fW.fit=T,init.val=init)
W.fit.34 = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,locs=par.locs,pen.const=20,method="BFGS",joint.fit=F,fW.fit=T,init.val=init/init[1])
# save.image(file="d3pollution_fits_tau95_bww_penconst_search.RData")


xstar.sims.cond = sim.cond(w=wexc,r0w=r0wexc,nsim=50000,k=1,
                           gfun=function(w,par,locs=par.locs){gfun.pwl(w,par,locs)},par=R.fit.24$mle)
xstar.sims.joint = sim.joint(nsim=50000,k.vals=1,shape=3,
                             par.locs=par.locs,par.locs.W=par.locs,
                             par=R.fit.24$mle, fW.par=W.fit.34$fW.mle,
                             bww=qr$bww, tau=tau, r=r, w=w)[[1]]

save.image(file="d3pollution_fit.RData")
stop()

##############################################################################

# # experiment with the L1 penalty
# R.fit.13.L1pen = fit.pwlin(r=rexc,r0w=r0wexc,w=wexc,locs=par.locs,
#                            pen.const.L1=0.1,pen.const.L2=0.0,method="BFGS",
#                            init.val=rep(0.5,dim(par.locs)[1]))#init.val=init)
# mle =  R.fit.13.L1pen$mle
# g.vals = gfun.3d(x=qr$wpts,par=mle, ref.angles=par.locs)
# plot3d(ds.exp.3d/log(nrow(ds.exp.3d)))
# axes3d()
# title3d(xlab=colnames(ds.exp.3d)[1],ylab=colnames(ds.exp.3d)[2],zlab=colnames(ds.exp.3d)[3])
# surface3d(x=qr$wpts[,1]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
#           y=qr$wpts[,2]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
#           z=qr$wpts[,3]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
#           col="blue",alpha=0.4)
# view3d(userMatrix = usermat,zoom=0.8)

##############################################

# DIAGNOSTICS....
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
library(geometricMVE)
library(scales)

library(PWLExtremes)

library(this.path)
setwd(this.path::here())
load("d3pollution_fit.RData")
# load("d3pollution_fits_tau95_betterlambda.RData")
# load("d3pollution_fits_finermesh.RData")

# pen.consts = sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$pen))

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

# load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")
# names(pen.fit.res[[1]])
# 
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$pp.score.13)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$pp.score.24)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$pp.score.5)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$pp.score.6)))
# 
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$return.score.13)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$return.score.24)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$return.score.5)))
# which.min(sapply(1:length(pen.fit.res), function(ii) return(pen.fit.res[[ii]]$return.score.6)))
# 
# nnodes=nrow(par.locs)

plot3d(ds.exp.3d,col="black",alpha=0.2)
par3d(windowRect=c(0,0, 725, 725))
surface3d(x=qr$wpts[,1]*matrix(qr$r.tau.wpts,30,30),
          y=qr$wpts[,2]*matrix(qr$r.tau.wpts,30,30),
          z=qr$wpts[,3]*matrix(qr$r.tau.wpts,30,30),
          col="red",alpha=0.5)
axes3d(edges="bbox")
view3d(userMatrix = usermat,zoom=0.8)
snapshot3d("~/Desktop/d3pollution_thresh.png",width=500,height=500)

plot3d(ds.exp.3d,col="black",alpha=0.2)
par3d(windowRect=c(0,0, 725, 725))
points3d(ds.exp.3d[excind,],col="blue",size=5)
axes3d(edges="bbox")
view3d(userMatrix = usermat,zoom=0.8)
snapshot3d("~/Desktop/d3pollution_exc.png",width=500,height=500)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_thresh.png",width=500,height=500)


# qr$wpts = cbind(qr$wpts,1-apply(qr$wpts,1,sum))
mle =  R.fit.24$mle#pen.fit.res[[1]]$RW.fit.5$mle
g.vals = gfun.pwl(x=qr$wpts,par=mle, ref.angles=par.locs)
plot3d(ds.exp.3d/log(nrow(ds.exp.3d)))
# k1 <- max(qr$wpts[, 1] * qr$r.tau.wpts, na.rm = TRUE)
# k2 <- max(qr$wpts[, 2] * qr$r.tau.wpts, na.rm = TRUE)
# k3 <- max(qr$wpts[, 3] * qr$r.tau.wpts, na.rm = TRUE)
# surface3d(x=qr$wpts[,1]*qr$r.tau.wpts/k1,
#           y=qr$wpts[,2]*qr$r.tau.wpts/k2,
#           z=qr$wpts[,3]*qr$r.tau.wpts/k3,
#           col="red",alpha=0.4)
axes3d()
title3d(xlab=colnames(ds.exp.3d)[1],ylab=colnames(ds.exp.3d)[2],zlab=colnames(ds.exp.3d)[3])
surface3d(x=qr$wpts[,1]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
          y=qr$wpts[,2]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
          z=qr$wpts[,3]/matrix(g.vals,nrow=sqrt(length(g.vals)),ncol=sqrt(length(g.vals))),
          col="blue",alpha=0.4)
view3d(userMatrix = usermat,zoom=0.8)
snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_g_est.png",width=500,height=500)

# check g(x,y,1)=1
obj = function(xy){ gfun.pwl(x=c(xy,1),par=mle, ref.angles=par.locs) }
optim.out = optim(par=rep(0.5,2),obj)
round(optim.out$par,3)
gfun.pwl(x=c(optim.out$par,1),par=mle, ref.angles=par.locs)
gfun.pwl(x=c(1,1,1),par=mle, ref.angles=par.locs)

# PP-QQ plots
gfun = function(w,par,locs=par.locs){
  gfun.pwl(w,par,locs)
}
num <- pgamma(rexc, shape = 3, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
den <- pgamma(r0wexc, shape = 3, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
nr <- length(rexc)
load("d3pollution_chi_BSfits.Rdata")
model.PPQQ.bs = do.call(rbind,lapply(bs.res.pwl, function(lst) {
    BS.mle = lst$pwl.bs.pars$R.fit.24.bs
    if(is.null(BS.mle)){
      return(NA)
    } else {
      num <- pgamma(rexc, shape = 3, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
      den <- pgamma(r0wexc, shape = 3, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
      return(sort(1 - num/den))
    }
}))
model.PP.CI.low = apply(model.PPQQ.bs,2,quantile,p=0.025,na.rm=T)
model.PP.CI.upp = apply(model.PPQQ.bs,2,quantile,p=0.975,na.rm=T)
model.QQ.CI.low = apply(qexp(model.PPQQ.bs),2,quantile,p=0.025,na.rm=T)
model.QQ.CI.upp = apply(qexp(model.PPQQ.bs),2,quantile,p=0.975,na.rm=T)

pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_PP.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL,xlab = "empirical", ylab = "model",
         xlim=c(0,1),ylim=c(0,1))
CI.x = c(c(1:nr)/(nr + 1),rev(c(1:nr)/(nr + 1)))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(model.PP.CI.low,rev(model.PP.CI.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0
polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(c(1:nr)/(nr + 1), sort(1 - num/den), pch = 20, cex = 0.8,
     xlim=c(0,1),ylim=c(0,1))
abline(a = 0, b = 1)
dev.off()
pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_QQ.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL,xlab = "empirical", ylab = "model", pch = 20, cex = 0.8,
     xlim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))),
     ylim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))))
CI.x = c(qexp(c(1:nr)/(nr + 1)),rev(qexp(c(1:nr)/(nr + 1))))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(model.QQ.CI.low,rev(model.QQ.CI.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0
polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(qexp(c(1:nr)/(nr + 1)), qexp(sort(1 - num/den)), pch = 20, cex = 0.8)
abline(a = 0, b = 1)
dev.off()

#########################################################

# chi plots

idx.groups = lapply(c(1:3), function(x) combn(3,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

thresh.r = as.numeric(qr$r.tau.wpts)
w.mesh = qr$wpts
u.vals.lst = lapply(idx.groups,function(idx){
  u.min = pexp(max(sapply(1:length(thresh.r),function(w.idx){
    w.vec = w.mesh[w.idx,]
    r0w.val = thresh.r[w.idx]
    return(r0w.val*min(w.vec[idx],na.rm=T))
  }),na.rm=T))
  return(seq(u.min,1,length.out=100))
})
# save(u.vals.lst,file="u_vals_d3pollution.RData")

# chi.vals.emp = list(list(NULL))
xstar.sims = xstar.sims.joint#sim.cond(w=wexc,r0w=r0wexc,nsim=50000,k=1,gfun=gfun,par=R.fit.24$mle)
col.names = as.character(c(1,2,3))
chi.vals = list(list(NULL))
grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = u.vals.lst[[grp.no]]
  
  # empirical
  chi.u.emp = sapply(u.vals, function(u) mean(apply(ds.exp.3d[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u))  # assuming the margins are approximately standard exponential
  
  # pwl
  chi.u.est.pwl = mean(excind,na.rm=T)*sapply(u.vals, function(u) mean(apply(xstar.sims[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
  
  col.names = as.character(c(1,2,3))
  col.names = col.names[which.cols]
  name = paste(col.names, collapse = "")
  
  chi.vals[[grp.no]] = list(u.vals=u.vals,
                            chi.u.emp=chi.u.emp,
                            chi.u.est.pwl=chi.u.est.pwl,
                            name=name)
  
  # ylab = bquote(.(rlang::sym("chi"))[.(name)](u))
  # 
  # pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_chi_plot_",name,".pdf"),height=5,width=5)
  # par(mfrow=c(1,1),pty="s",mar=c(4,4.5,1,4.5),cex.axis=1.5,cex.lab=1.5)
  # plot(u.vals,chi.u.emp,type="l",lwd=1,
  #      xlab="u",ylab=ylab,
  #      ylim=c(0,1.0))
  # lines(u.vals,chi.u.est.pwl,lty=2,col="blue",lwd=2)
  # dev.off()
  
  grp.no=grp.no+1
}

load("d3pollution_chi_bs_emp.Rdata")
grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = chi.vals[[grp.no]]$u.vals
  
  # empirical
  chi.u.emp = chi.vals[[grp.no]]$chi.u.emp
  emp.CI = do.call(rbind,lapply(bs.res.emp,function(lst){
    lst[[grp.no]]
  }))
  emp.CI.low = apply(emp.CI,2,quantile,p=0.025,na.rm=T)
  emp.CI.upp = apply(emp.CI,2,quantile,p=0.975,na.rm=T)
  
  # pwl
  chi.u.est.pwl = chi.vals[[grp.no]]$chi.u.est.pwl
  pwl.CI = do.call(rbind,lapply(bs.res.pwl,function(lst){
    lst$pwl.res$pwl.res4[[grp.no]]
  }))
  pwl.CI.low = apply(pwl.CI,2,quantile,p=0.025,na.rm=T)
  pwl.CI.upp = apply(pwl.CI,2,quantile,p=0.975,na.rm=T)
  
  col.names = as.character(c(1,3,4))
  col.names = col.names[which.cols]
  name = chi.vals[[grp.no]]$name
  
  ylab = bquote(.(rlang::sym("chi"))[.(name)](u))

  pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_chi_plot_",name,".pdf"),height=5,width=5)
  par(mfrow=c(1,1),pty="s",mar=c(4,4.5,1,4.5),cex.axis=1.5,cex.lab=1.5)
  plot(u.vals,chi.u.emp,type="l",lwd=1,
       xlab="u",ylab=ylab,
       ylim=c(0,1.0))
  CI.x = c(u.vals,rev(u.vals))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(emp.CI.low,rev(emp.CI.upp))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=alpha("gray",0.3),border=NA)
  lines(u.vals,chi.u.est.pwl,lty=2,col="blue",lwd=2)
  CI.x = c(u.vals,rev(u.vals))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(pwl.CI.low,rev(pwl.CI.upp))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=alpha("blue",0.25),border=NA)
  dev.off()
  
  grp.no=grp.no+1
}

# Return levels
gauge.est.vals.w = gfun.pwl(w,par=R.fit.24$mle,ref.angles=par.locs)
# T.vals = seq(ceil(1/(1-tau)),1000,length.out=100)+1
T.vals = seq(21,1000,by=10)
est.prop.exc = sapply(T.vals,function(T.val){
  return.radii = qgamma(1-(1/T.val), shape=3, rate=gauge.est.vals.w)
  return(mean(r>return.radii))
})
# true.prop.exc = 1/T.vals

# x.vals = true.prop.exc[-c(1:10)]
# y.vals = est.prop.exc[-c(1:10)]
x.vals = log(T.vals)
y.vals = log(1/est.prop.exc)

est.log.T.BS = do.call(rbind,lapply(bs.res.pwl,function(lst) {
  BS.mle = lst$pwl.bs.pars$R.fit.24.bs
  if(is.null(BS.mle)){
    return(NA)
  } else {
    BS.gauge.est.vals.w = gfun.pwl(w,par=BS.mle,ref.angles=par.locs)
    return(log(1/sapply(T.vals,function(T.val){
      return.radii = qgamma(1-(1/T.val), shape=3, rate=BS.gauge.est.vals.w)
      return(mean(r>return.radii))
    })))
  }
}))
est.log.T.BS.CI.low = apply(est.log.T.BS,2,quantile,p=0.025,na.rm=T)
est.log.T.BS.CI.upp = apply(est.log.T.BS,2,quantile,p=0.975,na.rm=T)

pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_return.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL,
     xlim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     ylim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     xlab="log(T), true",ylab="log(T), estimated",type="p")
CI.x = c(x.vals,rev(x.vals))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(est.log.T.BS.CI.low,rev(est.log.T.BS.CI.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0
polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(x=x.vals,y=y.vals,pch=20)
segments(-1,-1,10000,10000)
# lines(x.vals,pwl.CI.l[-c(1:10)],lty=2)
# lines(x.vals,pwl.CI.u[-c(1:10)],lty=2)
dev.off()

# plot a few return level set boundaries
T.vals=c(50,100,1000)
gauge.est.vals = gfun.pwl(qr$wpts,par=R.fit.24$mle, ref.angles=par.locs)
for(T.val in T.vals){
  return.radii = qgamma(1-(1/T.val)*(1/(1-tau)), shape=3, rate=gauge.est.vals)
  return.radii = matrix(return.radii,nrow=sqrt(nrow(qr$wpts)))
  
  open3d()
  points3d(ds.exp.3d)
  axes3d(edges="bbox")
  rgl::surface3d(qr$wpts[, 1] *return.radii,
                 qr$wpts[, 2] *return.radii,
                 qr$wpts[, 3] *return.radii,
                 col = "grey", alpha = 0.4, add = T)
  title3d(xlab=colnames(ds.exp.3d)[1],
          ylab=colnames(ds.exp.3d)[2],
          zlab=colnames(ds.exp.3d)[3])
  view3d(userMatrix = usermat,zoom=0.8)
  snapshot3d(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d3pollution_return_T",T.val,".png"),width=500,height=500)
}

