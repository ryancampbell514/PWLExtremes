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

## @knitr data
# load in the data (exponential margins)
load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")
x.E = ds.exp.3d
x.F = qfrechet(pexp(x.E))

load("d3pollution_fit.RData")

# plot3d(x.E)
# plot3d(x.F)

# Interested in the probability that PM10 is large when others are small,
# this probability *should* be positive
C = 3
notC = c(1:3)[-C]
u.val.exp = KDE.quant.eval(wpts = matrix(c(0,0,1),nrow=1),r=r,w=w,ker="Gaussian",bww=0.075)
u = qfrechet(pexp(u.val.exp))
q.seq = seq(u+0.01,qfrechet(0.999),length.out=50)  # a sequence above u
delta = 0.3

# plot the region in exponential margins
plot3d(x.E,pch=20)
xyz.mesh = expand.grid(replicate(3, seq(0, max(x.E), len = 30), simplify = FALSE))
cond = xyz.mesh[,C] > qexp(pfrechet(u)) & apply(xyz.mesh[,notC],1,max) <= xyz.mesh[,C]^delta
points3d(xyz.mesh[cond,],col="red")
view3d(userMatrix = usermat)


## @knitr probest
# Method 1: empirical
prob.full.emp = sapply(q.seq, function(qq){
  mean(x.E[,C]>(-log(1-exp(-1/qq))) & apply(x.E[,notC],1,max)<=((x.E[,C])^delta))
})

# Method 2: Simpson et al.
prob.decond = mean(apply(x.F[,notC],1,max)<=(x.F[,C])^delta)  
q.vals = x.F[,C]
tau.est = min(sum((q.vals>u)*log(q.vals/u))/sum(q.vals>u),1)
K.C = mean(q.vals>u) * u^(1/tau.est)
prob.cond = K.C*(q.seq^(-1/tau.est))
prob.full = prob.cond*prob.decond

# Method 3: geometric
load("d3pollution_fit.RData")
x.star = xstar.sims.cond  # here, we perform *slightly* better on conditional sample (no angular density)
prob.full.geom = mean(excind,na.rm=T) *
  sapply(q.seq, function(qq){
    mean(x.star[,C]>(-log(1-exp(-1/qq))) & apply(x.star[,notC],1,max)<=((x.star[,C])^delta))
  })

# plot the estimates
par(pty="s",mfrow=c(1,1))
plot(qexp(pfrechet(q.seq)),prob.full.emp,type="l",lwd=2,xlab="q",ylim=c(0,max(prob.full.emp,prob.full,prob.full.geom)))
lines(qexp(pfrechet(q.seq)),prob.full.geom,col="blue",lwd=2,lty=2)
lines(qexp(pfrechet(q.seq)),prob.full,col="red",lwd=2,lty=3)
legend("topright",lwd=rep(2,3),col=c("black","blue","red"),legend=c("emp","PWL","Simp."),
       lty=c(1,2,3))
abline(h=0)
