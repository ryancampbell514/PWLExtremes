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
# load("d4pollution_fit_morenodes.RData")
load("../d4pollution_fit_akdethresh.RData")
w.full=w
r.full=r

source("../../../../../geometricMVE_RCcode/diagnostics_new.R")

plotfittedthresh.3dproj(qr,which.proj = 1)
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj1_akde.html")
webshot2::webshot("threshproj1_akde.html", "threshproj1_akde.png", vwidth = 500, vheight = 500)
# file.remove("fittedgauge_proj_1.html")

plotfittedthresh.3dproj(qr,which.proj = 2)
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj2_akde.html")
webshot2::webshot("threshproj2_akde.html", "threshproj2_akde.png", vwidth = 500, vheight = 500)

plotfittedthresh.3dproj(qr,which.proj = 3)
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj3_akde.html")
webshot2::webshot("threshproj3_akde.html", "threshproj3_akde.png", vwidth = 500, vheight = 500)

plotfittedthresh.3dproj(qr,which.proj = 4)
view3d(userMatrix = usermat,zoom=0.8)
htmlwidgets::saveWidget(rglwidget(),file="threshproj4_akde.html")
webshot2::webshot("threshproj4_akde.html", "threshproj4_akde.png", vwidth = 500, vheight = 500)