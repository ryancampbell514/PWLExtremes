rm(list = ls())

library(geometry)
library(pracma)
library(openair)
library(lubridate)
library(lattice)
# library(geometricMVE)
library(scales)

Sys.setenv(RGL_USE_NULL = TRUE)
library(rgl)
options(rgl.printRglwidget = TRUE)

library(this.path)
setwd(this.path::here())

fn.dir = "../../../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# library(PWLExtremes)

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

# load("d4pollution_fit.RData")  # <-- this doesn't exist
# load("d4pollution_chi_BSfits.Rdata")
load("../d4pollution_fit_morenodes.RData")
w.full=w
r.full=r

# source("../../../../../geometricMVE_RCcode/diagnostics_new.R")
# source("../../../../../geometricMVE_RCcode/thresh_estimation.R")
source("../../../../../geometricMVEdev/thresh_estimation_AKDE.R")
source("../../../../../geometricMVEdev/diagnostics_new.R")

qr = fit.thresh(r=r,w=w,tau=tau,bww=0.075)

# modify old qr to match new syntax
# qr$w = qr$w[,-4]
# qr$method="KDE"
# qr$r = r.full
# qr$bww = matrix(qr$bww,nrow=nrow(qr$w),ncol=ncol(qr$w))

nms = colnames(w)

plotfittedthresh.3dproj(qr,which.proj = 1, xlab = nms[2], ylab = nms[3], zlab = nms[4])
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj1_orig.html")
webshot2::webshot("threshproj1_orig.html", "threshproj1_orig.png", vwidth = 500, vheight = 500)
# file.remove("fittedgauge_proj_1.html")

plotfittedthresh.3dproj(qr,which.proj = 2, xlab = nms[1], ylab = nms[3], zlab = nms[4])
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj2_orig.html")
webshot2::webshot("threshproj2_orig.html", "threshproj2_orig.png", vwidth = 500, vheight = 500)

plotfittedthresh.3dproj(qr,which.proj = 3, xlab = nms[1], ylab = nms[2], zlab = nms[4])
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj3_orig.html")
webshot2::webshot("threshproj3_orig.html", "threshproj3_orig.png", vwidth = 500, vheight = 500)

plotfittedthresh.3dproj(qr,which.proj = 4, xlab = nms[1], ylab = nms[2], zlab = nms[3])
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj4_orig.html")
webshot2::webshot("threshproj4_orig.html", "threshproj4_orig.png", vwidth = 500, vheight = 500)