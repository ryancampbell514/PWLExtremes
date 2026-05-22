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

qr.akde = fit.thresh(r=r,w=w,tau=tau,bww=0.075)
qr = qr.akde
r0w = qr$r0w
excind<-r>r0w
print(mean(excind,na.rm=T))

for(a in seq(0,1,by=0.25)){
  qr.akde = fit.thresh(r=r,w=w,tau=tau,alpha=a)
  qr = qr.akde
  r0w = qr$r0w
  excind<-r>r0w
  print(mean(excind,na.rm=T))
}

load("../fitting_orig/d4pollution_fit_morenodes.RData")
mean(excind)
