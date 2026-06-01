rm(list=ls())

library(rgl)
library(geometricMVE)

############################################################

# d=2

w.mesh = seq(0,1,length.out=1000)

gfun1 = function(xy,par=0.5){
  max((xy[1]-xy[2])/par, (xy[2]-xy[1])/par, sum(xy)/(2-par))
}
gvals1 = apply(cbind(w.mesh,1-w.mesh,0),1,gfun1)

par(pty="s",mar=c(2,2,1,2))
plot(cbind(w.mesh,1-w.mesh)/gvals1,type="l",lwd=2,xlim=c(0,1.2),ylim=c(0,1.2),
     xlab=" ",ylab=" ",axes = F)
axis(side=1)
axis(side=2)
segments(1,0,1,1,lty=3)
segments(0,1,1,1,lty=3)

gfun2 = function(xy){
  geometricMVE::gauge_logistic(xy,par=0.4)
}
gvals2 = apply(cbind(w.mesh,1-w.mesh),1,gfun2)

par(pty="s",mar=c(2,2,1,2))
plot(cbind(w.mesh,1-w.mesh)/gvals2,type="l",lwd=2,xlim=c(0,1.2),ylim=c(0,1.2),
     xlab=" ",ylab=" ",axes = F)
axis(side=1)
axis(side=2)
segments(1,0,1,1,lty=3)
segments(0,1,1,1,lty=3)

############################################################

# d=3

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T)  #par3d()$userMatrix

num.angles=500
w.mesh.fine = seq(0,1,length.out=num.angles)
w.mesh.fine = expand.grid(w.mesh.fine,w.mesh.fine)
w.mesh.fine = cbind(w.mesh.fine,1-apply(w.mesh.fine,1,sum))
w.mesh.fine[,3] = ifelse(w.mesh.fine[,3]<0,NA,w.mesh.fine[,3])
names(w.mesh.fine)=c("w1","w2","w3")

gfun1 = function(xyz){
  geometricMVE::gauge_alogistic(vec=xyz,par=rep(0.4,3),
                                 theta1=0,theta2=0,theta3=0,theta12=1,theta13=1,theta23=1,theta123=0)
}
g.vals.fine = apply(w.mesh.fine,1,gfun1)
g.vals.fine.mat = matrix(unlist(g.vals.fine),num.angles,num.angles)

open3d()
surface3d(w.mesh.fine[,1]/g.vals.fine.mat,
          w.mesh.fine[,2]/g.vals.fine.mat,
          w.mesh.fine[,3]/g.vals.fine.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]))
view3d(userMatrix = usermat,zoom=0.8)

gfun2 = function(xyz){
  geometricMVE::gauge_alogistic(vec=xyz,par=rep(0.4,1),
                                 theta1=0,theta2=0,theta3=0,theta12=0,theta13=0,theta23=0,theta123=1)
}
g.vals.fine = apply(w.mesh.fine,1,gfun2)
g.vals.fine.mat = matrix(unlist(g.vals.fine),num.angles,num.angles)

open3d()
surface3d(w.mesh.fine[,1]/g.vals.fine.mat,
          w.mesh.fine[,2]/g.vals.fine.mat,
          w.mesh.fine[,3]/g.vals.fine.mat,
          col="blue",alpha=0.5)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]))
view3d(userMatrix = usermat,zoom=0.8)
