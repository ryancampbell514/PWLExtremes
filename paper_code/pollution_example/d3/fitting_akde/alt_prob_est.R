## @knitr setup

# d=3 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
library(geometricMVE)
library(scales)

library(this.path)
setwd(this.path::here())

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
library(PWLExtremes)

library(evd)

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

############################################################################

# load in the data (exponential margins)
load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")

# load in the fitted models
load("d3pollution_fit.RData")

# Interested in the probability that PM10 is large when others are small,
# Want to capture the tail decau of this
quantile.surface.cart = qr$wpts * qr$r.tau.wpts
which.max.pm10.loc = which.max(quantile.surface.cart[,3])
u0 = pexp(quantile.surface.cart[which.max.pm10.loc,3])
u.seq = seq(u0,1,length.out=100)  # a sequence above u

delta1 = pexp(8)
delta2 = pexp(2)

## @knitr probest
# Method 1: empirical
prob.full.emp = sapply(u.seq, function(u.val){
  mean(ds.exp.3d[,1]<qexp(delta1) & ds.exp.3d[,2]<qexp(delta2) &
         ds.exp.3d[,3]>qexp(u.val))
})

bs.est.res.emp = NULL
num.bs.reps = 500
for(rep in c(1:num.bs.reps)){
  set.seed(1000+rep)
  
  message(paste("bootstrap rep:",rep))
  
  # block bootstrap
  # acf(ds.exp.3d)
  block.size = 7
  n.idx.sample = nrow(ds.exp.3d)
  idx.init = sample(c(1:nrow(ds.exp.3d)),ceil(nrow(ds.exp.3d)/block.size),replace=T)
  idx = lapply(idx.init,function(ii) ii+c(0:(block.size-1)))
  idx = unlist(idx)
  idx = idx[idx<=nrow(ds.exp.3d)]
  # idx = sample(c(1:nrow(ds.exp.3d)),nrow(ds.exp.3d),replace=T)   # regular bootstrap
  ds.bs = ds.exp.3d[idx,]
  # acf(ds.bs)
  
  emp.res = sapply(u.seq, function(u.val){
    mean(ds.bs[,1]<qexp(delta1) & ds.bs[,2]<qexp(delta2) &
           ds.bs[,3]>qexp(u.val))
  })
  bs.est.res.emp[[rep]] = emp.res
}
bs.res.emp.mat = do.call(rbind,bs.est.res.emp)
emp.CI.low = apply(bs.res.emp.mat,2,quantile,0.025)
emp.CI.upp = apply(bs.res.emp.mat,2,quantile,0.975)

# Method 2: conditional extremes
X.exp=ds.exp.3d
source("~/Dropbox/phd_research/geomMVE_work/r_code/Functions/CEproflik.r")
th<-quantile(X.exp[,3],0.95)
## 1|3
fce1<-optim(HT.nll.PL,par=c(0.5,0.5),dat=X.exp[,c(1,3)],u=th,which=2)
Z1<-(X.exp[X.exp[,3]>th,1]-fce1$par[1]*X.exp[X.exp[,3]>th,3])/X.exp[X.exp[,3]>th,3]^fce1$par[2]
## 2|3
fce2<-optim(HT.nll.PL,par=c(0.5,0.5),dat=X.exp[,c(2,3)],u=th,which=2)
Z2<-(X.exp[X.exp[,3]>th,2]-fce2$par[1]*X.exp[X.exp[,3]>th,3])/X.exp[X.exp[,3]>th,3]^fce2$par[2]

Z<-cbind(Z1,Z2)
n1<-50000
xe.star<-rexp(n1)+qexp(u0)
star.ind<-sample(1:length(Z1),size=n1,replace=T)
Zstar<-Z[star.ind,]
Y1.star<-fce1$par[1]*xe.star+Zstar[,1]*xe.star^fce1$par[2]
Y2.star<-fce2$par[1]*xe.star+Zstar[,2]*xe.star^fce2$par[2]
xstar.CE<-cbind(Y1.star,Y2.star,xe.star)
# plot3d(X.exp)
# points3d(xstar.CE,color="red")
prob.full.CE = exp(-qexp(u0)) *sapply(u.seq, function(u.val){
  mean(xstar.CE[,1]<qexp(delta1) & xstar.CE[,2]<qexp(delta2) &
         xstar.CE[,3]>qexp(u.val))
})

# Method 3: geometric
load("d3pollution_fit.RData")
x.star = xstar.sims.cond  # here, we perform *slightly* better on conditional sample (no angular density)
prob.full.geom = mean(excind,na.rm=T) *
  sapply(u.seq, function(u.val){
    mean(x.star[,1]<qexp(delta1) & x.star[,2]<qexp(delta2) &
           x.star[,3]>qexp(u.val))
  })

load("d3pollution_chi_BSfits.Rdata")
bs.est.res.pwl = NULL
num.bs.reps = 500
for(rep in c(1:num.bs.reps)){
  set.seed(1000+rep)
  
  block.size = 7
  n.idx.sample = nrow(ds.exp.3d)
  idx.init = sample(c(1:nrow(ds.exp.3d)),ceil(nrow(ds.exp.3d)/block.size),replace=T)
  idx = lapply(idx.init,function(ii) ii+c(0:(block.size-1)))
  idx = unlist(idx)
  idx = idx[idx<=nrow(ds.exp.3d)]
  
  ds.bs = ds.exp.3d[idx,]
  r.bs<-apply(ds.bs,1,sum)
  w.bs<-ds.bs/r.bs
  r0w.bs = qr$r0w[idx]
  excind.bs = excind[idx]
  rexc.bs<-r.bs[excind.bs]
  wexc.bs<-w.bs[excind.bs,]
  r0wexc.bs<-r0w.bs[excind.bs]
  na.ind.bs<-which(is.na(rexc.bs))
  # na.ind
  if(length(na.ind.bs)>0){
    rexc.bs<-rexc.bs[-na.ind.bs]
    wexc.bs<-wexc.bs[-na.ind.bs,]
    r0wexc.bs<-r0wexc.bs[-na.ind.bs]}
  
  message(paste("bootstrap rep:",rep))
  
  mle.bs = bs.res.pwl[[rep]]$pwl.bs.pars$R.fit.24.bs
  if(is.null(mle.bs)){
    next
  }
  xstar.sims.cond.bs = sim.cond(w=wexc.bs,r0w=r0wexc.bs,nsim=50000,k=1,
                                gfun=function(w,par,locs=par.locs){gfun.pwl(w,par,locs)},
                                par=mle.bs)
  
  pwl.res = mean(excind,na.rm=T) *
    sapply(u.seq, function(u.val){
      mean(xstar.sims.cond.bs[,1]<qexp(delta1) & xstar.sims.cond.bs[,2]<qexp(delta2) &
             xstar.sims.cond.bs[,3]>qexp(u.val))
    })
  bs.est.res.pwl[[rep]] = pwl.res
}
bs.res.pwl.mat = do.call(rbind,bs.est.res.pwl)
pwl.CI.low = apply(bs.res.pwl.mat,2,quantile,0.025,na.rm=T)
pwl.CI.upp = apply(bs.res.pwl.mat,2,quantile,0.975,na.rm=T)

# plot the estimates
pdf("~/Dropbox/phd_research/Lancs_Nice_seminar/figures/alt_prob_est_plot.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s",mar=c(4.2,4.2,4.2,4.2),cex.lab=1.5)
# par(pty="s",mfrow=c(1,1))
plot(u.seq,prob.full.emp,type="l",lwd=2,xlab="u",
     ylab=NA,
     ylim=c(0,max(prob.full.emp,prob.full.geom,emp.CI.upp,pwl.CI.upp)))
CI.x = c(u.seq,rev(u.seq))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(emp.CI.low,rev(emp.CI.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0
polygon(x=CI.x,y=CI.y,col=alpha("gray",0.3),border=NA)
lines(u.seq,prob.full.geom,col="blue",lwd=2,lty=2)
CI.x = c(u.seq,rev(u.seq))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(pwl.CI.low,rev(pwl.CI.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0
polygon(x=CI.x,y=CI.y,col=alpha("blue",0.25),border=NA)
abline(h=0)
# lines(u.seq,prob.full.CE,col="red",lwd=2,lty=3)
dev.off()
