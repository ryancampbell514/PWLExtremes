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

# fn.dir = "../../../../geometricMVE/R"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
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

############################################################################
mle = R.fit$mle


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

# est.log.T.BS = do.call(rbind,lapply(bs.res.pwl,function(lst) lst$est.log.T.bs))
# est.log.T.BS.CI.low = apply(est.log.T.BS,2,quantile,p=0.025)
# est.log.T.BS.CI.upp = apply(est.log.T.BS,2,quantile,p=0.975)

pdf("d4pollution_return.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL,
     xlim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     ylim=c(min(x.vals,y.vals),max(x.vals,y.vals)),
     xlab="log(T), true",ylab="log(T), estimated",type="p")
# CI.x = c(x.vals,rev(x.vals))
# CI.x[length(CI.x)+1] = CI.x[1]
# CI.y = c(est.log.T.BS.CI.low,rev(est.log.T.BS.CI.upp))
# CI.y[length(CI.y)+1] = CI.y[1]
# CI.y[is.na(CI.y)] = 0
# polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(x=x.vals,y=y.vals,pch=20)
segments(-1,-1,10000,10000)
# lines(x.vals,pwl.CI.l[-c(1:10)],lty=2)
# lines(x.vals,pwl.CI.u[-c(1:10)],lty=2)
dev.off()

