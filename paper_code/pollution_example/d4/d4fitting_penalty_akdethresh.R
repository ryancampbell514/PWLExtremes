# d=4 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
# library(geometricMVE)

# library(this.path)
# setwd(this.path::here())
setwd("/home/rcampbell/lancaster/pw_lin_gauge/pollution_example/d4")

# ##########################
# 
# # library(PWLExtremes)
# fn.dir = "../../../geometricMVE/R"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# source("../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")
# 
# ###########################
# 
# load("../ds_4d_urban.Rdata")
# head(ds.exp.4d)
# 
# r<-apply(ds.exp.4d,1,sum)
# w<-ds.exp.4d/r
# # #
# tau=0.9#0.70
# qr.akde = fit.thresh(r=r,w=w,tau=tau)
# qr = qr.akde
# r0w = qr$r0w
# excind<-r>r0w
# rexc<-r[excind]
# wexc<-w[excind,]
# r0w<-r0w[excind]
# save.image(file="d4pollution_fit_akdethresh.RData")
# 
#
# # load("d4pollution_fit.RData")
# par.locs = seq(0,1,by=0.2)  #0.2 or 1/6
# par.locs = as.matrix(expand.grid(par.locs,par.locs,par.locs))
# par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
# par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
# par.locs[,4] = ifelse(par.locs[,4]<0.001,0,par.locs[,4])
# par.locs = rbind(par.locs,
#                  rep(1/4,4),
#                 c(1/3,1/3,1/3,0),c(1/3,1/3,0,1/3),c(1/3,0,1/3,1/3),c(0,1/3,1/3,1/3),
#                 c(0.5,0.5,0,0),c(0.5,0,0.5,0),c(0.5,0,0,0.5),c(0,0.5,0.5,0),c(0,0.5,0,0.5),c(0,0,0.5,0.5))
# par.locs = as.matrix(par.locs[!duplicated(par.locs), ])
# par.locs = rbind(diag(4),
#                  rep(1/4,4),
#                  c(1/3,1/3,1/3,0),c(1/3,1/3,0,1/3),c(1/3,0,1/3,1/3),c(0,1/3,1/3,1/3),
#                  c(0.5,0.5,0,0),c(0.5,0,0.5,0),c(0.5,0,0,0.5),c(0,0.5,0.5,0),c(0,0.5,0,0.5),c(0,0,0.5,0.5))
# tri.mat = geometry::delaunayn(p=par.locs[,-4], output.options=TRUE)$tri
# new.mesh = t(apply(tri.mat,1,function(row.idx){
#   return(apply(par.locs[row.idx,],2,mean))
# }))
# new.mesh = new.mesh + rnorm(prod(dim(new.mesh)),sd=0.001)
# new.mesh = seq(0,1,by=0.2)  #seq(0,1,by=0.2)
# new.mesh = as.matrix(expand.grid(replicate(4-1, new.mesh, simplify = FALSE)))
# new.mesh = cbind(new.mesh,1-apply(new.mesh,1,sum))
# new.mesh = new.mesh / apply(new.mesh,1,sum)
# new.mesh = new.mesh[apply(new.mesh,1,function(w) !any(w<0)),]
# new.mesh = new.mesh[apply(new.mesh,1,function(vec) sum(vec>0)>2),]
# par.locs.2 = data.frame(rbind(par.locs,new.mesh))
# par.locs.2 = as.matrix(par.locs.2[!duplicated(par.locs.2), ])
# save(par.locs.2,file="parlocs2.RData")
# load("parlocs2.RData")
# par.locs=par.locs.2
# # nnodes=nrow(par.locs)
# #
# # # init = KDE.quant.eval(wpts=par.locs,r=r,w=w,tau=tau,bww=qr$bww,ker="Gaussian")
# init=eval.thresh(qr,par.locs)
# init = init/max(par.locs[,1] * init, na.rm = TRUE)
# # # init=NULL
# #
# message("Fitting R|W bounded")
# R.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
#                              fixshape=T,pen.const=1,
#                              method="BFGS",bound.fit=TRUE)
# print(R.fit$mle)
# save.image(file="d4pollution_fit_akdethresh.RData")
# 
# message("Fitting W")
# W.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
#                              fixshape=T,pen.const=1, #0.137931034482759,
#                              method="BFGS",W.fit=TRUE,joint.fit=FALSE)
# save.image(file="d4pollution_fit_akdethresh.RData")

load("d4pollution_fit_akdethresh.RData")

message("Fitting (R,W)")
RW.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
                          fixshape=T,pen.const=1,
                          method="BFGS",W.fit=TRUE,joint.fit=TRUE,bound.fit=TRUE)
print(RW.fit$mle)

# message("Simulating from models.")
# xstar.sims.cond = PWLExtremes::sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,
#                                         gfun=function(w,par,locs=par.locs){PWLExtremes::gfun.pwl(w,par,locs)},
#                                         par=R.fit$mle)
# save.image(file="d4pollution_fit_akdethresh.RData")

# TO DO: update sim.joint
source("../../sim_joint.R")
# xstar.sims.joint1 = sim.joint(nsim=50000,k.vals=1,shape=4,
#                               par.locs=par.locs,par.locs.W=par.locs,
#                               par=R.fit$mle, fW.par=W.fit$mle.W,
#                               fitted.thresh=qr)[[1]]
# save.image(file="d4pollution_fit_akdethresh.RData")
xstar.sims.joint2 = sim.joint(nsim=50000,k.vals=1,shape=4,
                             par.locs=par.locs,par.locs.W=par.locs,
                             par=RW.fit$mle, fW.par=RW.fit$mle.W,
                             fitted.thresh=qr)[[1]]
save.image(file="d4pollution_fit_akdethresh.RData")
