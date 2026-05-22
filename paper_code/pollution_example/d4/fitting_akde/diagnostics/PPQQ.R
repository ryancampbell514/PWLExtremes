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

#######################################################################


# PP-QQ plots
mle = R.fit$mle

gfun = function(w,par,locs=par.locs){
  gfun.pwl(w,par,locs)
}
num <- pgamma(rexc, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
den <- pgamma(r0w, shape = 4, rate = apply(wexc, 1, gfun, par = mle), lower.tail = F)
nr <- length(rexc)

# model.PPQQ.bs = do.call(rbind,lapply(bs.res.pwl, function(lst) {
#   BS.mle = lst$pwl.bs.pars$R.fit.24.bs
#   num <- pgamma(rexc, shape = 4, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
#   den <- pgamma(r0w, shape = 4, rate = apply(wexc, 1, gfun, par = BS.mle), lower.tail = F)
#   return(sort(1 - num/den))
# }))
# model.PP.CI.low = apply(model.PPQQ.bs,2,quantile,p=0.025)
# model.PP.CI.upp = apply(model.PPQQ.bs,2,quantile,p=0.975)
# model.QQ.CI.low = apply(qexp(model.PPQQ.bs),2,quantile,p=0.025)
# model.QQ.CI.upp = apply(qexp(model.PPQQ.bs),2,quantile,p=0.975)

pdf("d4pollution_akdethresh_PP.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL, xlab = "empirical", ylab = "model", xlim=c(0,1),ylim=c(0,1))
# CI.x = c(c(1:nr)/(nr + 1),rev(c(1:nr)/(nr + 1)))
# CI.x[length(CI.x)+1] = CI.x[1]
# CI.y = c(model.PP.CI.low,rev(model.PP.CI.upp))
# CI.y[length(CI.y)+1] = CI.y[1]
# CI.y[is.na(CI.y)] = 0
# polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(c(1:nr)/(nr + 1), sort(1 - num/den), pch = 20, cex = 0.8)
abline(a = 0, b = 1)
dev.off()
pdf("d4pollution_akdethresh_QQ.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
plot(NULL, xlab = "empirical", ylab = "model",
     xlim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))),
     ylim=c(min(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den))), max(qexp(c(1:nr)/(nr + 1)),qexp(sort(1 - num/den)))))
# CI.x = c(qexp(c(1:nr)/(nr + 1)),rev(qexp(c(1:nr)/(nr + 1))))
# CI.x[length(CI.x)+1] = CI.x[1]
# CI.y = c(model.QQ.CI.low,rev(model.QQ.CI.upp))
# CI.y[length(CI.y)+1] = CI.y[1]
# CI.y[is.na(CI.y)] = 0
# polygon(x=CI.x,y=CI.y,col=alpha("gray",0.8),border=alpha("gray",0.8))
points(qexp(c(1:nr)/(nr + 1)), qexp(sort(1 - num/den)),  pch = 20, cex = 0.8)
abline(a = 0, b = 1)
dev.off()