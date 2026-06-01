rm(list=ls())
library(geometricMVE)

source("extra-functions.R")

n.mesh = 500
w.mesh = seq(0,1,length.out=n.mesh)

### gaussian example
nnodes.vec = c(5,7,9)
gfun.true = apply(cbind(w.mesh,1-w.mesh),1,geometricMVE::gauge_gaussian,par=0.8)
set.seed(11)
x<-qexp(pnorm(rmvnorm(5000,sigma=matrix(c(1,0.8,0.8,1),2,2))))
for(nn in nnodes.vec){
  ref.angles = seq(0,1,length.out=nn)
  par.vals = 1/apply(cbind(ref.angles,1-ref.angles),1,geometricMVE::gauge_gaussian,par=0.8)
  # par.vals = par.vals/max(ref.angles*par.vals)
  gfun.pwl = apply(cbind(w.mesh,1-w.mesh),1,gfun.2d,par=par.vals,ref.angles=ref.angles)
  
  par(pty="s",mar=c(2,2,1,2))
  plot(x/log(5000),xlim=c(0,1.1),ylim=c(0,1.1),pch=20,cex=0.8,xlab=NA,ylab=NA,
       col=scales::alpha("grey",0.7))
  segments(1,0,1,1,lty=3)
  segments(0,1,1,1,lty=3)
  lines(cbind(w.mesh,1-w.mesh)/gfun.true,lwd=2)
  lines(cbind(w.mesh,1-w.mesh)/gfun.pwl,lwd=2,col="blue",lty=2)
}