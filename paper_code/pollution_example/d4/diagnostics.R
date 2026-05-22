rm(list = ls())

library(geometry)
library(pracma)
library(openair)
library(lubridate)
library(lattice)
# library(geometricMVE)
library(scales)

library(rgl)
options(rgl.printRglwidget = TRUE)

library(this.path)
setwd(this.path::here())

fn.dir = "../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# library(PWLExtremes)

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

# load("d4pollution_fit.RData")  # <-- this doesn't exist
# load("d4pollution_chi_BSfits.Rdata")
# load("d4pollution_fit_morenodes.RData")
load("d4pollution_fit_akdethresh.RData")
w.full=w
r.full=r

###############################################################################

# Plot the threshold

# Function to project from 4 to 3d
# tau=0.70
# thr=fit.thresh(r=r.full,w=w.full,tau=tau,bww=0.075)
source("../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")
thr=fit.thresh(r=r.full,w=w.full,tau=tau)

rtilde.proj<-function(x,which.min,upper=10,res=100)
{
  if(any(is.na(x))){
    return(NA)
  }
  # x is a single angle
  if(which.min==1){
    dummy<-function(y)
    {
      w.inp = c(y,x)
      sum(w.inp) / eval.thresh(thr,w.eval = matrix(w.inp/sum(w.inp),nrow=1))
    }
  } else if(which.min==2){
    dummy<-function(y)
    {
      w.inp = c(x[1],y,x[2:3])
      sum(w.inp) / eval.thresh(thr,w.eval = matrix(w.inp/sum(w.inp),nrow=1))
    }
  } else if(which.min==3){
    dummy<-function(y)
    {
      w.inp = c(x[1:2],y,x[3])
      sum(w.inp) / eval.thresh(thr,w.eval = matrix(w.inp/sum(w.inp),nrow=1))
    }
  }else if(which.min==4){
    dummy<-function(y)
    {
      w.inp = c(x,y)
      sum(w.inp) / eval.thresh(thr,w.eval = matrix(w.inp/sum(w.inp),nrow=1))
    }
  }
  opt<-optimize(dummy,interval=c(0,upper))
  return(opt$objective)
}

n.mesh = 30
wseq<-seq(0,1,len=n.mesh)
wgrid = expand.grid(wseq,wseq)
wmat<-cbind(wgrid,1-apply(wgrid,1,sum))
nacond = wmat[,3] <0

thresh.proj1<-apply(wmat,1,rtilde.proj,which.min=1)
thresh.proj1[nacond] = NA
thresh.proj1.mat = matrix(thresh.proj1,sqrt(length(thresh.proj1)),sqrt(length(thresh.proj1)))
open3d()
surface3d(wmat[,1]/thresh.proj1.mat,
          wmat[,2]/thresh.proj1.mat,
          wmat[,3]/thresh.proj1.mat,
          col="red",alpha=0.4)
points3d(ds.exp.4d[,-1],col="black",alpha=0.3)
axes3d()
view3d(userMatrix = usermat)
htmlwidgets::saveWidget(rglwidget(),file="threshproj1_akde.html")
snapshot3d("threshproj1.png",width=500,height=500)

thresh.proj2<-apply(wmat,1,rtilde.proj,which.min=2)
thresh.proj2[nacond] = NA
thresh.proj2.mat = matrix(thresh.proj2,sqrt(length(thresh.proj2)),sqrt(length(thresh.proj2)))
open3d()
surface3d(wmat[,1]/thresh.proj2.mat,
          wmat[,2]/thresh.proj2.mat,
          wmat[,3]/thresh.proj2.mat,
          col="red",alpha=0.4)
points3d(ds.exp.4d[,-2],col="black",alpha=0.3)
axes3d()
view3d(userMatrix = usermat)
htmlwidgets::saveWidget(rglwidget(),file="threshproj2_akde.html")
# snapshot3d("threshproj2.png",width=500,height=500)

thresh.proj3<-apply(wmat,1,rtilde.proj,which.min=3)
thresh.proj3[nacond] = NA
thresh.proj3.mat = matrix(thresh.proj3,sqrt(length(thresh.proj3)),sqrt(length(thresh.proj3)))
open3d()
surface3d(wmat[,1]/thresh.proj3.mat,
          wmat[,2]/thresh.proj3.mat,
          wmat[,3]/thresh.proj3.mat,
          col="red",alpha=0.4)
points3d(ds.exp.4d[,-3],col="black",alpha=0.3)
axes3d()
view3d(userMatrix = usermat)
htmlwidgets::saveWidget(rglwidget(),file="threshproj3_akde.html")
# snapshot3d("threshproj3.png",width=500,height=500)

thresh.proj4<-apply(wmat,1,rtilde.proj,which.min=4)
thresh.proj4[nacond] = NA
thresh.proj4.mat = matrix(thresh.proj4,sqrt(length(thresh.proj4)),sqrt(length(thresh.proj4)))
open3d()
surface3d(wmat[,1]/thresh.proj4.mat,
          wmat[,2]/thresh.proj4.mat,
          wmat[,3]/thresh.proj4.mat,
          col="red",alpha=0.4)
points3d(ds.exp.4d[,-4],col="black",alpha=0.3)
axes3d()
view3d(userMatrix = usermat)
htmlwidgets::saveWidget(rglwidget(),file="threshproj4_akde.html")
# snapshot3d("threshproj4.png",width=500,height=500)

save.image("akde_thresh_proj.RData")

stop()

###############################################################################

n.mesh = 100
mle = R.fit$mle

# print the componentwise max
w.mesh <- expand.grid(replicate(4 - 1, seq(0, 1, len = n.mesh), simplify = FALSE))
w.mesh = cbind(w.mesh, 1-apply(w.mesh,1,sum))
w.mesh = w.mesh[!apply(w.mesh,1,function(vec) any(vec<0)),]
gvals = gfun.pwl(w.mesh,mle,par.locs)
print(apply(w.mesh / gvals, 2, max, na.rm=T))

# What's the minimum of g(x1,x2,1,x4)?
obj = function(xyz){ gfun.pwl(x=c(xyz[1:2],1,xyz[3]),par=mle, ref.angles=par.locs) }
optim.out = optim(par=rep(0.5,3),obj)
round(optim.out$par,3)
optim.out$value
# gfun.pwl(x=c(optim.out$par,1),par=mle, ref.angles=par.locs)

## METHOD 1: minimise over a mesh
proj.g.fn = function(gfun,w,which.w,...){
  # gfun -> gauge that takes in 4-dim vectors
  # w -> 3-min input
  # which.w -> which index to take min over
  
  w.inp = matrix(NA,nrow=n.mesh,ncol=4)
  w.inp[,-which.w] = matrix(as.numeric(w),ncol=3,nrow=n.mesh,byrow=T)
  w.inp[,which.w] = seq(0,1,length.out=n.mesh)
  
  # return(min(apply(w.inp,1,gfun,...)))
  return(min(gfun(w.inp,...)))
}

wpts<-expand.grid(seq(0,1,len=n.mesh),seq(0,1,len=n.mesh))
wpts = cbind(wpts,1-apply(wpts,1,sum))

t1 = Sys.time()
gvals.est1 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=1,par=mle,ref.angles=par.locs)
gvals.est1.mat = matrix(gvals.est1,n.mesh,n.mesh)
t2 = Sys.time()
print(t2-t1)

gvals.est2 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=2,par=mle,ref.angles=par.locs)
gvals.est2.mat = matrix(gvals.est2,n.mesh,n.mesh)

gvals.est3 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=3,par=mle,ref.angles=par.locs)
gvals.est3.mat = matrix(gvals.est3,n.mesh,n.mesh)

gvals.est4 = apply(wpts,1,proj.g.fn,gfun=gfun.pwl,which.w=4,par=mle,ref.angles=par.locs)
gvals.est4.mat = matrix(gvals.est4,n.mesh,n.mesh)

save.image(file="d4pollution_gauge_projections_AKDEthresh.RData")
# stop()


denom = log(nrow(ds.exp.4d))
ds.scaled = ds.exp.4d / denom

open3d()
plot3d(ds.scaled[,c(2,3,4)])
surface3d(wpts[,1]/gvals.est1.mat,
          wpts[,2]/gvals.est1.mat,
          wpts[,3]/gvals.est1.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_1.png",width=500,height=500)
# htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_1.html")

open3d()
plot3d(ds.scaled[,c(1,3,4)])
surface3d(wpts[,1]/gvals.est2.mat,
          wpts[,2]/gvals.est2.mat,
          wpts[,3]/gvals.est2.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_2.png",width=500,height=500)
# htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_2.html")

open3d()
plot3d(ds.scaled[,c(1,2,4)])
surface3d(wpts[,1]/gvals.est3.mat,
          wpts[,2]/gvals.est3.mat,
          wpts[,3]/gvals.est3.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_3.png",width=500,height=500)
# htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_3.html")

open3d()
plot3d(ds.scaled[,c(1,2,3)])
surface3d(wpts[,1]/gvals.est4.mat,
          wpts[,2]/gvals.est4.mat,
          wpts[,3]/gvals.est4.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_4.png",width=500,height=500)
# htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_4.html")

###############################################################################

# PP-QQ plots
mle = R.fit$mle

gfun = function(w,par,locs=par.locs){
  gfun.pwl(w,par,locs)
}
num <- pgamma(rexc, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
den <- pgamma(r0w, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
nr <- length(rexc)

# model.PPQQ.bs = do.call(rbind,lapply(bs.res.pwl, function(lst) {
#   BS.mle = lst$pwl.bs.pars$R.fit.24.bs
#   num <- pgamma(rexc, shape = 4, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
#   den <- pgamma(r0w, shape = 4, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
#   return(sort(1 - num/den))
# }))
# model.PP.CI.low = apply(model.PPQQ.bs,2,quantile,p=0.025)
# model.PP.CI.upp = apply(model.PPQQ.bs,2,quantile,p=0.975)
# model.QQ.CI.low = apply(qexp(model.PPQQ.bs),2,quantile,p=0.025)
# model.QQ.CI.upp = apply(qexp(model.PPQQ.bs),2,quantile,p=0.975)

pdf("d4pollution_akdethresh_PP.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL, xlab = "empirical", ylab = "model", xlim=c(0,1),ylim=c(0,1))
# CI.x = c(c(1:nr)/(nr + 1),rev(c(1:nr)/(nr + 1)))
# CI.x[length(CI.x)+1] = CI.x[1]
# CI.y = c(model.PP.CI.low,rev(model.PP.CI.upp))
# CI.y[length(CI.y)+1] = CI.y[1]
# CI.y[is.na(CI.y)] = 0
# polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(c(1:nr)/(nr + 1), sort(1 - num/den), pch = 20, cex = 0.8)
abline(a = 0, b = 1)
dev.off()
pdf("d4pollution_akdethresh_QQ.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL, xlab = "empirical", ylab = "model",
     xlim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))),
     ylim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))))
# CI.x = c(qexp(c(1:nr)/(nr + 1)),rev(qexp(c(1:nr)/(nr + 1))))
# CI.x[length(CI.x)+1] = CI.x[1]
# CI.y = c(model.QQ.CI.low,rev(model.QQ.CI.upp))
# CI.y[length(CI.y)+1] = CI.y[1]
# CI.y[is.na(CI.y)] = 0
# polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(qexp(c(1:nr)/(nr + 1)), qexp(sort(1 - num/den)),  pch = 20, cex = 0.8)
abline(a = 0, b = 1)
dev.off()

#############################################################

# chi plots

idx.groups = lapply(c(1:4), function(x) combn(4,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

thresh.r = as.numeric(qr$r.tau.wpts)
w.mesh = qr$wpts
# u.vals.lst = lapply(idx.groups,function(idx){
#   u.min = pexp(max(sapply(1:length(thresh.r),function(w.idx){
#     w.vec = w.mesh[w.idx,]
#     r0w.val = thresh.r[w.idx]
#     return(r0w.val*min(w.vec[idx],na.rm=T))
#   }),na.rm=T))
#   return(seq(u.min,1,length.out=100))
# })
# save(u.vals.lst,file="u_vals_d4pollution.RData")
load("u_vals_d4pollution.RData")

# chi.vals.emp = list(list(NULL))
# gfun = function(w,par,locs=par.locs){
#   gfun.pwl(w,par,locs)
# }
xstar.sims = xstar.sims.joint#xstar.sims.joint#sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,gfun=gfun,par=mle)
col.names = as.character(c(1:4))
chi.vals = list(list(NULL))
grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = u.vals.lst[[grp.no]]
  
  # empirical
  chi.u.emp = sapply(u.vals, function(u) mean(apply(ds.exp.4d[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u))  # assuming the margins are approximately standard exponential
  
  # pwl
  chi.u.est.pwl = mean(excind,na.rm=T)*sapply(u.vals, function(u) mean(apply(xstar.sims[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
  
  col.names = as.character(c(1:4))
  col.names = col.names[which.cols]
  name = paste(col.names, collapse = "")
  
  chi.vals[[grp.no]] = list(u.vals=u.vals,
                            chi.u.emp=chi.u.emp,
                            chi.u.est.pwl=chi.u.est.pwl,
                            name=name)
  
  # ylab = bquote(.(rlang::sym("chi"))[.(name)](u))
  # 
  # # pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_chi_plot_",name,".pdf"),height=5,width=5)
  # par(mfrow=c(1,1),pty="s",mar=c(4,4.5,1,4.5),cex.axis=1.5,cex.lab=1.5)
  # plot(u.vals,chi.u.emp,type="l",lwd=1,
  #      xlab="u",ylab=ylab,
  #      ylim=c(0,1.0))
  # lines(u.vals,chi.u.est.pwl,lty=2,col="blue",lwd=2)
  # # dev.off()
  
  grp.no=grp.no+1
}

grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = chi.vals[[grp.no]]$u.vals
  chi.u.emp = chi.vals[[grp.no]]$chi.u.emp
  emp.CI = do.call(rbind,lapply(bs.res.pwl,function(lst){
    lst$pwl.res$emp.res[[grp.no]]
  }))
  emp.CI.low = apply(emp.CI,2,quantile,p=0.025,na.rm=T)
  emp.CI.upp = apply(emp.CI,2,quantile,p=0.975,na.rm=T)
  chi.u.est.pwl = chi.vals[[grp.no]]$chi.u.est.pwl
  pwl.CI =  do.call(rbind,lapply(bs.res.pwl,function(lst){
    lst$pwl.res$pwl.res[[grp.no]]
  }))
  pwl.CI.low = apply(pwl.CI,2,quantile,p=0.025,na.rm=T)
  pwl.CI.upp = apply(pwl.CI,2,quantile,p=0.975,na.rm=T)
    
  col.names = as.character(c(1:4))
  col.names = col.names[which.cols]
  name = chi.vals[[grp.no]]$name
  
  ylab = bquote(.(rlang::sym("chi"))[.(name)](u))

  pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_chi_plot_",name,".pdf"),height=5,width=5)
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


#############################################################

# Return level score

gauge.est.vals.w = gfun.pwl(w,par=mle,ref.angles=par.locs)
# T.vals = seq(ceil(1/(1-tau)),1000,length.out=100)+1
T.vals = seq(10,1000,by=10)
est.prop.exc = sapply(T.vals,function(T.val){
  return.radii = qgamma(1-(1/T.val), shape=4, rate=gauge.est.vals.w)
  return(mean(r>return.radii))
})
# true.prop.exc = 1/T.vals

# x.vals = true.prop.exc[-c(1:10)]
# y.vals = est.prop.exc[-c(1:10)]
x.vals = log(T.vals)
y.vals = log(1/est.prop.exc)

est.log.T.BS = do.call(rbind,lapply(bs.res.pwl,function(lst) lst$est.log.T.bs))
est.log.T.BS.CI.low = apply(est.log.T.BS,2,quantile,p=0.025)
est.log.T.BS.CI.upp = apply(est.log.T.BS,2,quantile,p=0.975)

pdf("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_return.pdf",width=5,height=5)
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

