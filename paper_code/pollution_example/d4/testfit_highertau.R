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

# library(PWLExtremes)
fn.dir = "../../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
source("../../AKDE_new_offsimplex_threshold/thresh_estimation_constrained_AKDE.R")

###########################

load("../ds_4d_urban.Rdata")
head(ds.exp.4d)

r<-apply(ds.exp.4d,1,sum)
w<-ds.exp.4d/r
# # 
tau=0.9
qr.akde = fit.thresh(r=r,w=w,tau=tau)
qr = qr.akde
r0w = qr$r0w
excind<-r>r0w
rexc<-r[excind]
wexc<-w[excind,]
r0w<-r0w[excind]

load("parlocs2.RData")
par.locs=par.locs.2

init=eval.thresh(qr,par.locs)
init = init/max(par.locs[,1] * init, na.rm = TRUE)

message("Fitting R|W bounded")
R.fit = fit.geometric.pwl(r=rexc,r0w=r0w,w=wexc,locs=par.locs,init.val=init,
                          fixshape=T,pen.const=1, 
                          method="BFGS",bound.fit=TRUE)

mle = R.fit$mle

gfun = function(w,par,locs=par.locs){
  gfun.pwl(w,par,locs)
}
num <- pgamma(rexc, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
den <- pgamma(r0w, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
nr <- length(rexc)

par(mfrow=c(1,2),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL, xlab = "empirical", ylab = "model", xlim=c(0,1),ylim=c(0,1))
points(c(1:nr)/(nr + 1), sort(1 - num/den), pch = 20, cex = 0.8)
abline(a = 0, b = 1)
plot(NULL, xlab = "empirical", ylab = "model",
     xlim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))),
     ylim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))))
points(qexp(c(1:nr)/(nr + 1)), qexp(sort(1 - num/den)),  pch = 20, cex = 0.8)
abline(a = 0, b = 1)