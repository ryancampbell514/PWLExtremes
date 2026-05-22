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

options(rgl.printRglwidget = TRUE)

##########################

# load-in my version of geometricMVE
fn.dir = "../../../../geometricMVE_RCcode/"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# load-in some needed scripts that aren't on geometricMVE
source("../../../../geometricMVEdev/thresh_estimation_AKDE.R")
source("../../../../geometricMVEdev/sim_joint.R")

###########################

# Load in the data (exponential margins)
load("../../ds_4d_urban.Rdata")
head(ds.exp.4d)
print(dim(ds.exp.4d))

# convert to radii and angles
r<-apply(ds.exp.4d,1,sum)
w<-ds.exp.4d/r

# fit the threshold
tau=0.70

# hw.mat = ks::hpi(w[,-c(2:4)])
# ?ks::kde

# bww = apply(w[,-4],2,function(ww){
#   pilot <- stats::density(ww, bw = "nrd0")
#   fhat_i <- stats::approx(pilot$x, pilot$y, xout = ww)$y
#   g <- exp(mean(log(fhat_i)))   # geometric mean
#   lambda <- sqrt(g / fhat_i)
#   h_global <- stats::bw.nrd0(ww)
#   h_i <- h_global * lambda
# })

qr.akde = fit.thresh(r=r,w=w,tau=tau,bww=NULL,alpha=0)# 0.75 gave bad chi plots
# plotfittedthresh.3dproj(qr.akde, resolution=30, add=FALSE, which.proj=1)

# obtain the exceedances
qr = qr.akde
r0w = qr$r0w
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind,]
r0w<-r0w[excind]
save.image(file="d4pollution_fit_akdethresh.RData")

# load-in the reference angles
load("../parlocs2.RData")
par.locs=par.locs.2
rm(par.locs.2)
nnodes=nrow(par.locs)

# the initial values for optimisation
init = eval.thresh(qr,par.locs)
init = init/max(par.locs[,1] * init, na.rm = TRUE)

# Fit the models
message("fitting R")
R.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
                             fixshape=T,pen.const=1,bound.fit=TRUE)
xstar.sims.cond = PWLExtremes::sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,
                                        gfun=function(w,par,locs=par.locs){PWLExtremes::gfun.pwl(w,par,locs)},
                                        par=R.fit$mle)
save.image(file="d4pollution_fit_akdethresh.RData")

message("Fitting W")
W.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
                             fixshape=T,pen.const=1, #0.137931034482759,
                             method="BFGS",W.fit=TRUE,joint.fit=FALSE)
save.image(file="d4pollution_fit_akdethresh.RData")

message("Fitting (R,W)")
RW.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
                          fixshape=T,pen.const=1,
                          method="BFGS",W.fit=TRUE,joint.fit=TRUE,bound.fit=TRUE)
save.image(file="d4pollution_fit_akdethresh.RData")

###################################

# draw simulations
message("Simulating from models.")
# xstar.sims.cond = PWLExtremes::sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,
#                                         gfun=function(w,par,locs=par.locs){PWLExtremes::gfun.pwl(w,par,locs)},
#                                         par=R.fit$mle)
# save.image(file="d4pollution_fit_akdethresh.RData")

xstar.sims.joint2 = sim.joint(nsim=50000,k.vals=1,shape=4,
                             par.locs=par.locs,par.locs.W=par.locs,
                             par=RW.fit$mle, fW.par=RW.fit$mle.W,
                             fitted.thresh=qr)[[1]]
save.image(file="d4pollution_fit_akdethresh.RData")
