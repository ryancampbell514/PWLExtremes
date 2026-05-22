# d=4 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
# library(geometricMVE)

library(this.path)
setwd(this.path::here())

# library(PWLExtremes)
fn.dir = "../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
source("../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")


# # ###########################
# # 
# load("../ds_4d_urban.Rdata")
# head(ds.exp.4d)
# # 
# # # # check AD and AI groups
# # # source("~/Dropbox/phd_research/pw_lin_gauge/SimStudy_3d/par/simstudypar/Extra Functions/SimpsonEtAlFunctions.R")
# # # m31<-method2(ds.exp.4d, delta=0.4, hillQuantile=0.9, pi= 0.05)
# # # m31
# # # m32<-method2(ds.exp.4d, delta=0.5, hillQuantile=0.9, pi= 0.05)
# # # m32
# # # m33<-method2(ds.exp.4d, delta=0.6, hillQuantile=0.9, pi= 0.05)
# # # m33
# # # 
# # # thetavec1<-rep(NA,length(m31$prob))
# # # thetavec1[m31$prob>0]<-1
# # # thetavec1[m31$prob==0]<-10e10
# # # 
# # # thetavec2<-rep(NA,length(m32$prob))
# # # thetavec2[m32$prob>0]<-1
# # # thetavec2[m32$prob==0]<-10e10
# # # 
# # # thetavec3<-rep(NA,length(m32$prob))
# # # thetavec3[m33$prob>0]<-1
# # # thetavec3[m33$prob==0]<-10e10
# # 
# r<-apply(ds.exp.4d,1,sum)
# w<-ds.exp.4d/r
# # 
# tau=0.70  #0.95
# # # qr=radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=T,bww=0.05) #NULL
# # #qr=PWLExtremes::radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=FALSE,bww=0.0222222222222222,ker="Epanechnikov")
# # qr=PWLExtremes::radial.quants.L1.KDE(r=r,w=w,tau=tau,mesh.eval=TRUE,bww=0.075,ker="Gaussian")
# # # qr.res = list(ds=ds.exp.4d,
# # #               tau=tau,
# # #               qr=qr)
# # #save.image(file="d4pollution_fits_tau70_bww_penconst_search.RData")
# # # load("d4pollution_quants_tau070_bwsearch.RData")
# # #tau = qr.res$tau
# # #qr = qr.res$qr
# # # load("d4pollution_fits_tau70_bww_penconst_search.RData")
# qr.akde = fit.thresh(r=r,w=w,tau=tau)
# qr = qr.akde
# r0w = qr$r0w
# excind<-r>r0w
# rexc<-r[excind]
# wexc<-w[excind,]
# r0w<-r0w[excind]
# 
# # load("d4pollution_fit.RData")
# # par.locs = seq(0,1,by=1/6)  #0.2
# # par.locs = as.matrix(expand.grid(par.locs,par.locs,par.locs))
# # par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
# # par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
# # par.locs[,4] = ifelse(par.locs[,4]<0.001,0,par.locs[,4])
# # par.locs = rbind(diag(4),
# #                  rep(1/4,4),
# #                  c(1/3,1/3,1/3,0),c(1/3,1/3,0,1/3),c(1/3,0,1/3,1/3),c(0,1/3,1/3,1/3),
# #                  c(0.5,0.5,0,0),c(0.5,0,0.5,0),c(0.5,0,0,0.5),c(0,0.5,0.5,0),c(0,0.5,0,0.5),c(0,0,0.5,0.5))
# load("parlocs2.RData")
# par.locs=par.locs.2
# nnodes=nrow(par.locs)
# 
# # init = KDE.quant.eval(wpts=par.locs,r=r,w=w,tau=tau,bww=qr$bww,ker="Gaussian")
# init=eval.thresh(qr,par.locs)
# init = init/max(par.locs[,1] * init, na.rm = TRUE)
# # init=NULL
# 
# message("Fitting R|W bounded")
# # t1 = Sys.time()
# # R.fit.24 = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
# #                      pen.const=1, #0.137931034482759,
# #                      method="BFGS",bound.fit=TRUE)
# R.fit.24 = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
#                              fixshape=T,pen.const=1, #0.137931034482759,
#                              method="BFGS",bound.fit=TRUE)
# # t2 = Sys.time()
# # print(t2-t1)
# # save.image(file="d4pollution_fit_morenodes.RData")
# save.image(file="d4pollution_fit_akdethresh.RData")
# 
# message("Fitting W")
# t1 = Sys.time()
# # W.fit = fit.pwlin(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init/init[1],
# #                   pen.const=20,method="BFGS",fW.fit=TRUE)
# W.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
#                              fixshape=T,pen.const=1, #0.137931034482759,
#                              method="BFGS",W.fit=TRUE)
# t2 = Sys.time()
# print(t2-t1)
# # save.image(file="d4pollution_fit_morenodes.RData")
# save.image(file="d4pollution_fit_akdethresh.RData")

load("d4pollution_fit_akdethresh.RData")


xstar.sims.cond = PWLExtremes::sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,
                           gfun=function(w,par,locs=par.locs){PWLExtremes::gfun.pwl(w,par,locs)},
                           par=R.fit.24$mle)
xstar.sims.joint = PWLExtremes::sim.joint(nsim=50000,k.vals=1,shape=4,
                             par.locs=par.locs,par.locs.W=par.locs,
                             par=R.fit.24$mle, fW.par=W.fit$fW.mle,
                             bww=qr$bww, tau=tau, r=r, w=w)[[1]]

# save.image(file="d4pollution_fit_morenodes.RData")
save.image(file="d4pollution_fit_akdethresh.RData")

stop()

# Gauge function projections





# 
# 
# ##########################################
# 
# # re-fit at these MLEs,
# 
# # load("d4pollution_fits_tau90_bestlambda.RData")
# 
# tri.mat = geometry::delaunayn(p=par.locs[,-4], output.options=TRUE)$tri
# new.mesh = t(apply(tri.mat,1,function(row.idx){
#   return(apply(par.locs[row.idx,],2,mean))
# }))
# new.mesh = new.mesh + rnorm(prod(dim(new.mesh)),sd=0.001)
# # new.mesh = seq(0,1,by=0.2)  #seq(0,1,by=0.2)
# # new.mesh = as.matrix(expand.grid(replicate(4-1, new.mesh, simplify = FALSE)))
# # new.mesh = cbind(new.mesh,1-apply(new.mesh,1,sum))
# new.mesh = new.mesh / apply(new.mesh,1,sum)
# new.mesh = new.mesh[apply(new.mesh,1,function(w) !any(w<0)),]
# new.mesh = new.mesh[apply(new.mesh,1,function(vec) sum(vec>0)>2),]
# par.locs.2 = data.frame(rbind(par.locs,new.mesh))
# par.locs.2 = as.matrix(par.locs.2[!duplicated(par.locs.2), ])
# save(par.locs.2,file="parlocs2.RData")
# # load("parlocs2.RData")
# 
# init.2.mod13 = 1/gfun.pwl(x = par.locs.2, par = R.fit.13$mle, ref.angles = par.locs)
# # init.2.mod13 = 1/gfun.pwl(x = par.locs.2, par = res$R.fit.13$mle, ref.angles = par.locs)
# # init.2.mod24 = 1/gfun.pwl(x = par.locs.2, par = res$R.fit.24$mle, ref.angles = par.locs)
# # init.2.mod5 = 1/gfun.pwl(x = par.locs.2, par = res$RW.fit.5$mle, ref.angles = par.locs)
# # init.2.mod6 = 1/gfun.pwl(x = par.locs.2, par = res$RW.fit.6$mle, ref.angles = par.locs)
# 
# message("Re-fitting R|W unbounded")
# t1 = Sys.time()
# R.fit.13.refit = fit.pwlin.4d(r=rexc,r0w=r0w,w=wexc,locs=par.locs.2,init.val=init.2.mod13,pen.const=0.0001,#0.0005,#0.0002343750,
#                               method="BFGS")
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# # R.fit.13.refit$par.locs = par.locs.2
# # save(R.fit.13.refit,file="R_fit_13_refit.RData")
# 
# message("Re-fitting R|W bounded")
# t1 = Sys.time()
# R.fit.24.refit = bound.fit.4d(r=rexc,r0w=r0w,w=wexc,locs=par.locs.2,init.val=init.2.mod24,pen.const=0.0001,#0.003750000,
#                               method="BFGS")
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# message("Re-fitting R,W unbounded")
# t1 = Sys.time()
# RW.fit.5.refit = fit.pwlin.4d(r=rexc,r0w=r0w,w=wexc,locs=par.locs.2,init.val=init.2.mod5,pen.const=0.00005,#0.001023102,
#                               method="BFGS",joint.fit=T,fW.fit=T)
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# message("Re-fitting R,W bounded")
# t1 = Sys.time()
# RW.fit.6.refit = bound.fit.4d(r=rexc,r0w=r0w,w=wexc,locs=par.locs.2,init.val=init.2.mod6,pen.const=0.00005,#0.005424805,
#                               method="BFGS",joint.fit=T,fW.fit=T)
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# # W.fit.34 = fit.pwlin.4d(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,pen.const=pen.const*10,method="BFGS",joint.fit=F,fW.fit=T)
# 
# message("Sampling R|W unbounded")
# t1 = Sys.time()
# x.star.1.refit = sim.cond(w=wexc,r0w=r0w,nsim=50000,shape=4,par=R.fit.13.refit$mle,
#                     gfun=function(w,par,locs=par.locs.2){gfun.pwl(x=w,par=par,ref.angles=locs)})
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# message("Sampling R|W bounded")
# t1 = Sys.time()
# x.star.2.refit = sim.cond(w=wexc,r0w=r0w,nsim=50000,shape=4,par=R.fit.24.refit$mle,
#                     gfun=function(w,par,locs=par.locs.2){gfun.pwl(x=w,par=par,ref.angles=locs)})
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# message("Sampling R,W joint")
# t1 = Sys.time()
# x.star.3.refit = sim.joint(nsim=50000,k.vals=1,shape=4,par=R.fit.13.refit$mle,
#                            fW.par=W.fit$fW.mle,par.locs=par.locs.2,par.locs.W=par.locs,
#                            r=r,w=w,tau=0.7)
# t2 = Sys.time()
# print(t2-t1)
# 
# message("Sampling R,W joint, bounded")
# t1 = Sys.time()
# x.star.4.refit = sim.joint(nsim=50000,k.vals=1,shape=4,par=R.fit.24.refit$mle,
#                            fW.par=W.fit$fW.mle,par.locs=par.locs.2,par.locs.W=par.locs,
#                            r=r,w=w,tau=0.7)
# t2 = Sys.time()
# print(t2-t1)
# 
# message("Sampling R,W unbounded")
# t1 = Sys.time()
# x.star.5.refit = sim.joint(nsim=50000,
#                      # gfun=function(w,par,locs=par.locs.2){gfun.pwl(x=w,par=par,ref.angles=locs)},
#                      shape=4,par=RW.fit.5.refit$mle,fW.par=RW.fit.5.refit$mle,
#                      par.locs=par.locs.2,r=r,w=w,tau=tau)
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# message("Sampling R,W bounded")
# t1 = Sys.time()
# x.star.6.refit = sim.joint(nsim=50000,
#                      # gfun=function(w,par,locs=par.locs.2){gfun.pwl(x=w,par=par,ref.angles=locs)},
#                      shape=4,par=RW.fit.6.refit$mle,fW.par=RW.fit.6.refit$mle,
#                      par.locs=par.locs.2,r=r,w=w,tau=tau)
# t2 = Sys.time()
# print(t2-t1)
# save.image(file="d4pollution_fits_tau070_temp.RData")
# 
# res.refit = list(par.locs=par.locs,par.locs.2=par.locs.2,qr=qr,tau=tau,
#            R.fit.13=R.fit.13.refit,R.fit.24=R.fit.24.refit,
#            W.fit=W.fit,
#            RW.fit.5=RW.fit.5.refit,RW.fit.6=RW.fit.6.refit,
#            x.star.1=x.star.1.refit,
#            x.star.2=x.star.2.refit,
#            x.star.3=x.star.3.refit[[1]],
#            x.star.4=x.star.4.refit[[1]],
#            x.star.5=x.star.5.refit[[1]],
#            x.star.6=x.star.6.refit[[1]])
# save(res.refit,file="d4pollution_fits_tau070.RData")
# stop()
