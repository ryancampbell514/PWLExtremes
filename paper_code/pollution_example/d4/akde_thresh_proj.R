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
load("d4pollution_chi_BSfits.Rdata")
load("d4pollution_fit_morenodes.RData")
# load("d4pollution_fit_akdethresh.RData")
w.full=w
r.full=r

tau=0.70
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
