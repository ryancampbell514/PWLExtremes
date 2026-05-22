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

fn.dir = "../../../../geometricMVE/R"
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
load("../d4pollution_fit_akdethresh.RData")
w.full=w
r.full=r

################################################################


n.mesh = 100
mle = RW.fit$mle

# print the componentwise max
w.mesh <- expand.grid(replicate(4 - 1, seq(0, 1, len = n.mesh), simplify = FALSE))
w.mesh = cbind(w.mesh, 1-apply(w.mesh,1,sum))
w.mesh = w.mesh[!apply(w.mesh,1,function(vec) any(vec<0)),]
gvals = gfun.pwl(w.mesh,mle,par.locs)
# print(apply(w.mesh / gvals, 2, max, na.rm=T))

# What's the minimum of g(x1,x2,1,x4)?
# obj = function(xyz){ gfun.pwl(x=c(xyz[1:2],1,xyz[3]),par=mle, ref.angles=par.locs) }
# optim.out = optim(par=rep(0.5,3),obj)
# round(optim.out$par,3)
# optim.out$value
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

save.image(file="d4pollution_gauge_projections_AKDEthresh_RW.RData")
# stop()

##############################################################################

load("d4pollution_gauge_projections_AKDEthresh_RW.RData")

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
htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_1.html")

open3d()
plot3d(ds.scaled[,c(1,3,4)])
surface3d(wpts[,1]/gvals.est2.mat,
          wpts[,2]/gvals.est2.mat,
          wpts[,3]/gvals.est2.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_2.png",width=500,height=500)
htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_2.html")

open3d()
plot3d(ds.scaled[,c(1,2,4)])
surface3d(wpts[,1]/gvals.est3.mat,
          wpts[,2]/gvals.est3.mat,
          wpts[,3]/gvals.est3.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_3.png",width=500,height=500)
htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_3.html")

open3d()
plot3d(ds.scaled[,c(1,2,3)])
surface3d(wpts[,1]/gvals.est4.mat,
          wpts[,2]/gvals.est4.mat,
          wpts[,3]/gvals.est4.mat,col="blue",alpha=0.4)
axes3d()
view3d(userMatrix = usermat)
# snapshot3d("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_fittedgauge_proj_4.png",width=500,height=500)
htmlwidgets::saveWidget(rglwidget(),file="fittedgauge_proj_4.html")
