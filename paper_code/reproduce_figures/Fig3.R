rm(list=ls())
library(geometricMVE)

rho<-0.8
ww = seq(0,1,by=0.001)

par.locs = seq(0,1,length.out=5)
par.locs = cbind(par.locs,1-par.locs)
test.pars = 1/apply(par.locs,1,geometricMVE::gauge_gaussian,par=0.8)

par(pty="s",mar=c(4,4,1,4))
plot(NULL,xlab=expression(w[1]),ylab=expression(w[2]),
     xlim=c(0,1.05),ylim=c(0,1.05))
apply(par.locs,1,function(ww){
  segments(x0=0,y0=0,x1=ww[1]*100,y1=ww[2]*100,lwd=2.5,col="darkslategrey")
})

par(pty="s",mar=c(4,4,1,4))
plot(par.locs*test.pars,type="p",pch=20,cex=2,xlab=expression(w[1]),ylab=expression(w[2]),
     xlim=c(0,1.05),ylim=c(0,1.05))
arrows(x0=test.pars[2]*par.locs[2,1],y0=test.pars[2]*par.locs[2,2],
       x1=test.pars[1]*par.locs[1,1],y1=test.pars[1]*par.locs[1,2],lwd=4,col="darkgreen",lty=1)
arrows(x0=test.pars[3]*par.locs[3,1],y0=test.pars[3]*par.locs[3,2],
       x1=test.pars[2]*par.locs[2,1],y1=test.pars[2]*par.locs[2,2],lwd=4,col="darkgreen",lty=1)
arrows(x0=test.pars[4]*par.locs[4,1],y0=test.pars[4]*par.locs[4,2],
       x1=test.pars[3]*par.locs[3,1],y1=test.pars[3]*par.locs[3,2],lwd=4,col="darkgreen",lty=1)
arrows(x0=test.pars[5]*par.locs[5,1],y0=test.pars[5]*par.locs[5,2],
       x1=test.pars[4]*par.locs[4,1],y1=test.pars[4]*par.locs[4,2],lwd=4,col="darkgreen",lty=1)
text(x=0.05,y=0.32,labels=expression(theta["1"]~"w*"^"1"),cex=1.1,col="darkgreen")
text(x=0.37,y=0.96,labels=expression(theta["2"]~"w*"^"2"),cex=1.1,col="darkgreen")
text(x=0.9,y=0.95,labels=expression(theta["3"]~"w*"^"3"),cex=1.1,col="darkgreen")
text(x=0.95,y=0.25,labels=expression(theta["4"]~"w*"^"4"),cex=1.1,col="darkgreen")
text(x=0.32,y=0.05,labels=expression(theta["5"]~"w*"^"5"),cex=1.1,col="darkgreen")
points(par.locs*test.pars,pch=20,cex=2)

par(pty="s",mar=c(4,4,1,4))
plot(par.locs*test.pars,type="p",pch=20,cex=0,xlab=expression(w[1]),ylab=expression(w[2]),
     xlim=c(0,1.05),ylim=c(0,1.05))
segments(test.pars[1]*par.locs[1,1],
         test.pars[1]*par.locs[1,2],
         test.pars[2]*par.locs[2,1],
         test.pars[2]*par.locs[2,2],lwd=2)
segments(test.pars[2]*par.locs[2,1],
         test.pars[2]*par.locs[2,2],
         test.pars[3]*par.locs[3,1],
         test.pars[3]*par.locs[3,2],lwd=2)
segments(test.pars[3]*par.locs[3,1],
         test.pars[3]*par.locs[3,2],
         test.pars[4]*par.locs[4,1],
         test.pars[4]*par.locs[4,2],lwd=2)
segments(test.pars[4]*par.locs[4,1],
         test.pars[4]*par.locs[4,2],
         test.pars[5]*par.locs[5,1],
         test.pars[5]*par.locs[5,2],lwd=2)
polygon(rbind(test.pars*par.locs,0),col=scales::alpha("blue",0.4),lty=0)
segments(0,0,test.pars[1]*par.locs[1,1],test.pars[1]*par.locs[1,2],lty=2)
segments(0,0,test.pars[2]*par.locs[2,1],test.pars[2]*par.locs[2,2],lty=2)
segments(0,0,test.pars[3]*par.locs[3,1],test.pars[3]*par.locs[3,2],lty=2)
segments(0,0,test.pars[4]*par.locs[4,1],test.pars[4]*par.locs[4,2],lty=2)
segments(0,0,test.pars[5]*par.locs[5,1],test.pars[5]*par.locs[5,2],lty=2)
