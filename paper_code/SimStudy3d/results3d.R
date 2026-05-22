## @knitr results3d-setup

# Compile and plot d=2 simstudy results to show best model
rm(list=ls())
library(geometricMVE)
library(scales)

Sys.setenv(RGL_USE_NULL = TRUE)
library(rgl)
options(rgl.printRglwidget = TRUE)

library(geometry)

library(lattice)
library(squash)

# fn.dir = "/home/rcampbell/lancaster/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# setwd("/home/rcampbell/lancaster/pw_lin_gauge/SimStudy_3d/pwl_gradpen_finermesh_betterlambda/")
library(PWLExtremes)
library(this.path)
setwd(path.join(this.path::here(),"pwl_gitcode"))

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T)  #par3d()$userMatrix


##############################################################################

# ANGULAR DENSITY PLOTS

## Asymmetric Logistic 1
load("simstudy3/pwlinexp_SimStudy3d_alog1.Rdata")
wexc34.emp.samp = wexc
wexc34.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy5/pwlinexp_SimStudy3d_alog1.Rdata")
wexc5.emp.samp = wexc
wexc5.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy6/pwlinexp_SimStudy3d_alog1.Rdata")
wexc6.emp.samp = wexc
wexc6.mod.samp = xstar / apply(xstar,1,sum)

pdf("d3_alog1SS_Wexc_samp_marginal.pdf",width=9,height=12)
par(pty="s",mfrow=c(4,3),mar=c(4,4,4,4))
hist(wexc34.emp.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
dev.off()

## Asymmetric Logistic 2
load("simstudy3/pwlinexp_SimStudy3d_alog2.Rdata")
wexc34.emp.samp = wexc
wexc34.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy5/pwlinexp_SimStudy3d_alog2.Rdata")
wexc5.emp.samp = wexc
wexc5.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy6/pwlinexp_SimStudy3d_alog2.Rdata")
wexc6.emp.samp = wexc
wexc6.mod.samp = xstar / apply(xstar,1,sum)

pdf("d3_alog2SS_Wexc_samp_marginal.pdf",width=9,height=12)
par(pty="s",mfrow=c(4,3),mar=c(4,4,4,4))
hist(wexc34.emp.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
dev.off()


## Mixture Model
load("simstudy3/pwlinexp_SimStudy3d_mix.Rdata")
wexc34.emp.samp = wexc
wexc34.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy5/pwlinexp_SimStudy3d_mix.Rdata")
wexc5.emp.samp = wexc
wexc5.mod.samp = xstar / apply(xstar,1,sum)
load("simstudy6/pwlinexp_SimStudy3d_mix.Rdata")
wexc6.emp.samp = wexc
wexc6.mod.samp = xstar / apply(xstar,1,sum)

pdf("d3_mixSS_Wexc_samp_marginal.pdf",width=9,height=12)
par(pty="s",mfrow=c(4,3),mar=c(4,4,4,4))
hist(wexc34.emp.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.emp.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=20,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc34.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc5.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,1],freq = F,xlab=expression(w[1]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,2],freq = F,xlab=expression(w[2]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
hist(wexc6.mod.samp[,3],freq = F,xlab=expression(w[3]),main=NA,breaks=30,xlim=c(0,1),ylim=c(0,3.5))
dev.off()


##############################################################################

# GAUGE FUNCTION PLOTS

which.rep = c(1:200)

## @knitr results3d-gauges-alog1-ss1
load("simstudy1/SimStudy3d_alog1.Rdata")
pars.med = apply(pars[which.rep,],2,median,na.rm=T)
# w.mesh = qr$wpts
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
# w.mesh[w.mesh[,3]<0,]=rep(NA,3)
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh, byrow=F)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_pw3d(ww,par=rep(0.4,3))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)


open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog1_SS1.html")
# webshot2::webshot("alog1_SS1.html", "alog1_SS1.png", vwidth = 500, vheight = 500)
# file.remove("alog1_SS1.html")
snapshot3d("alog1_SS1.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog1-ss2
load("simstudy2/pwlinexp_SimStudy3d_alog1.Rdata")
pars.med = apply(pars[which.rep,],2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100) 
w.mesh = expand.grid(w.mesh,w.mesh) 
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_pw3d(ww,par=rep(0.4,3))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog1_SS2.html")
# webshot2::webshot("alog1_SS2.html", "alog1_SS2.png", vwidth = 500, vheight = 500)
# file.remove("alog1_SS2.html")
snapshot3d("alog1_SS2.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog1-ss2
load("simstudy5/pwlinexp_SimStudy3d_alog1.Rdata")
pars.med = apply(pars[which.rep,],2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_pw3d(ww,par=rep(0.4,3))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog1_SS5.html")
# webshot2::webshot("alog1_SS5.html", "alog1_SS5.png", vwidth = 500, vheight = 500)
# file.remove("alog1_SS5.html")
snapshot3d("alog1_SS5.png","png", width = 500, height = 600)


load("simstudy6/pwlinexp_SimStudy3d_alog1.Rdata")
pars.med = apply(pars[which.rep,],2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_pw3d(ww,par=rep(0.4,3))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog1_SS6.html")
# webshot2::webshot("alog1_SS6.html", "alog1_SS6.png", vwidth = 500, vheight = 500)
# file.remove("alog1_SS6.html")
snapshot3d("alog1_SS6.png","png", width = 500, height = 600)


## @knitr results3d-probs-alog1
rm(list=ls())

which.rep = c(1:200)

load("simstudy1/SimStudy3d_alog1.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3[which.rep]

load("simstudy2/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3[which.rep]

load("simstudy3/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss3 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3[which.rep]

load("simstudy4/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3[which.rep]

load("simstudy5/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3[which.rep]

load("simstudy6/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom[which.rep]
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2[which.rep]
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3[which.rep]
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom[which.rep]
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2[which.rep]
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3[which.rep]
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom[which.rep]
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2[which.rep]
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3[which.rep]

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_alog1.Rdata")

pdf("alog1_probests_main.pdf",width=7,height=3)
par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(prob1.est.geom2,pwl.prob1.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(prob2.est.geom2,pwl.prob2.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(prob3.est.geom3,pwl.prob3.est.geom3.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
dev.off()

pdf("alog1_probests_supp.pdf",width=9,height=9)
par(mfrow=c(4,1),pty="m",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(pwl.prob1.est.geom2.ss1,pwl.prob1.est.geom2.ss2,pwl.prob1.est.geom2.ss3,pwl.prob1.est.geom2.ss4,
        pwl.prob1.est.geom2.ss5,pwl.prob1.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(pwl.prob2.est.geom2.ss1,pwl.prob2.est.geom2.ss2,pwl.prob2.est.geom2.ss3,pwl.prob2.est.geom2.ss4,
        pwl.prob2.est.geom2.ss5,pwl.prob2.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
        pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
        pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[2,4]×[0.01,3] (no outliers)", outline=F)
abline(h=TP3,col="red")
dev.off()

########################################################################
########################################################################
rm(list=ls())

usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T)  #par3d()$userMatrix
usermat = par3d()$userMatrix

load("simstudy1/SimStudy3d_alog2.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_all3d(ww,par=rep(0.4,3),
                                   theta1=1,theta2=0,theta3=0,
                                   theta12=1,theta13=0,theta23=1,
                                   theta123=0)
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog2_SS1.html")
# webshot2::webshot("alog2_SS1.html", "alog2_SS1.png", vwidth = 500, vheight = 500)
# file.remove("alog2_SS1.html")
snapshot3d("alog2_SS1.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog2-ss2
load("simstudy2/pwlinexp_SimStudy3d_alog2.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_all3d(ww,par=rep(0.4,3),
                                   theta1=1,theta2=0,theta3=0,
                                   theta12=1,theta13=0,theta23=1,
                                   theta123=0)
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog2_SS2.html")
# webshot2::webshot("alog2_SS2.html", "alog2_SS2.png", vwidth = 500, vheight = 500)
# file.remove("alog2_SS2.html")
snapshot3d("alog2_SS2.png","png", width = 500, height = 600)

## @knitr results3d-gauges-alog2-ss2
load("simstudy5/pwlinexp_SimStudy3d_alog2.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_all3d(ww,par=rep(0.4,3),
                                   theta1=1,theta2=0,theta3=0,
                                   theta12=1,theta13=0,theta23=1,
                                   theta123=0)
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog2_SS5.html")
# webshot2::webshot("alog2_SS5.html", "alog2_SS5.png", vwidth = 500, vheight = 500)
# file.remove("alog2_SS5.html")
snapshot3d("alog2_SS5.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog2-ss2
load("simstudy6/pwlinexp_SimStudy3d_alog2.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    geometricMVE::gauge_rvad_all3d(ww,par=rep(0.4,3),
                                   theta1=1,theta2=0,theta3=0,
                                   theta12=1,theta13=0,theta23=1,
                                   theta123=0)
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="alog2_SS6.html")
# webshot2::webshot("alog2_SS6.html", "alog2_SS6.png", vwidth = 500, vheight = 500)
# file.remove("alog2_SS6.html")
snapshot3d("alog2_SS6.png","png", width = 500, height = 600)

rm(list=ls())
load("simstudy1/SimStudy3d_alog2.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3

load("simstudy2/pwlinexp_SimStudy3d_alog2.Rdata")
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3

# load("pwl_alog_N10/simstudy3/pwlinexp_SimStudy3d_alog2.Rdata")
load("simstudy3/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3

load("simstudy4/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3

load("simstudy5/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3

load("simstudy6/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_alog2.Rdata")

pdf("alog2_probests_main.pdf",width=7,height=3)
par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(prob1.est.geom2,pwl.prob1.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(prob2.est.geom2,pwl.prob2.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(prob3.est.geom3,pwl.prob3.est.geom3.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
dev.off()

pdf("alog2_probests_supp.pdf",width=9,height=7)
par(mfrow=c(3,1),pty="m",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(pwl.prob1.est.geom2.ss1,pwl.prob1.est.geom2.ss2,pwl.prob1.est.geom2.ss3,pwl.prob1.est.geom2.ss4,
        pwl.prob1.est.geom2.ss5,pwl.prob1.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(pwl.prob2.est.geom2.ss1,pwl.prob2.est.geom2.ss2,pwl.prob2.est.geom2.ss3,pwl.prob2.est.geom2.ss4,
        pwl.prob2.est.geom2.ss5,pwl.prob2.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
        pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
# boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
#         pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
#         names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
#         main="[8,10]×[2,4]×[0.01,3] (no outliers)", outline=F)
# abline(h=TP3,col="red")
dev.off()


##########################################################################
##########################################################################

rm(list=ls())
# source("/home/rcampbell/lancaster/pw_lin_gauge/Functions/gaugefunctions.R")
usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T)  #par3d()$userMatrix

load("simstudy1/SimStudy3d_mix.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    gauss.alog.mixgauge(ww,par=c(rep(0.4,3),rep(0.6,2)))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="mix_SS1.html")
# webshot2::webshot("mix_SS1.html", "mix_SS1.png", vwidth = 500, vheight = 500)
# file.remove("mix_SS1.html")
snapshot3d("mix_SS1.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog2-ss2
load("simstudy2/pwlinexp_SimStudy3d_mix.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    gauss.alog.mixgauge(ww,par=c(rep(0.4,3),rep(0.6,2)))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="mix_SS2.html")
# webshot2::webshot("mix_SS2.html", "mix_SS2.png", vwidth = 500, vheight = 500)
# file.remove("mix_SS2.html")
snapshot3d("mix_SS2.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog2-ss2
load("simstudy5/pwlinexp_SimStudy3d_mix.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    gauss.alog.mixgauge(ww,par=c(rep(0.4,3),rep(0.6,2)))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1.2)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="mix_SS5.html")
# webshot2::webshot("mix_SS5.html", "mix_SS5.png", vwidth = 500, vheight = 500)
# file.remove("mix_SS5.html")
snapshot3d("mix_SS5.png","png", width = 500, height = 600)


## @knitr results3d-gauges-alog2-ss2
load("simstudy6/pwlinexp_SimStudy3d_mix.Rdata")
pars.med = apply(pars,2,median,na.rm=T)
w.mesh = seq(0,1,length.out=100)
w.mesh = expand.grid(w.mesh,w.mesh)
w.mesh = cbind(w.mesh,1-apply(w.mesh,1,sum))
n.mesh = sqrt(nrow(w.mesh))
g.vals = apply(w.mesh,1,gfun,par=pars.med)
g.vals.mat = matrix(g.vals,n.mesh,n.mesh)
g.vals.true = apply(w.mesh,1,function(ww){
  if(any(ww<0)){
    return(NA)
  } else {
    gauss.alog.mixgauge(ww,par=c(rep(0.4,3),rep(0.6,2)))
  }
})
g.vals.true.mat = matrix(g.vals.true,n.mesh,n.mesh)

open3d()
rgl::material3d(depth_mask = FALSE)
surface3d(w.mesh[,1]/g.vals.mat,
          w.mesh[,2]/g.vals.mat,
          w.mesh[,3]/g.vals.mat,
          col="blue",alpha=0.4)
surface3d(w.mesh[,1]/g.vals.true.mat,
          w.mesh[,2]/g.vals.true.mat,
          w.mesh[,3]/g.vals.true.mat,
          col="red",alpha=0.4)
axes3d(edges="bbox")
title3d(xlab=expression(x[1]),ylab=expression(x[2]),zlab=expression(x[3]),cex=1)
view3d(userMatrix = usermat,zoom=0.8)
# htmlwidgets::saveWidget(rglwidget(),file="mix_SS6.html")
# webshot2::webshot("mix_SS6.html", "mix_SS6.png", vwidth = 500, vheight = 500)
# file.remove("mix_SS6.html")
snapshot3d("mix_SS6.png","png", width = 500, height = 600)


rm(list=ls())
load("simstudy1/SimStudy3d_mix.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3

load("simstudy2/pwlinexp_SimStudy3d_mix.Rdata")
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3

load("simstudy3/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3

load("simstudy4/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3

load("simstudy5/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3

load("simstudy6/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_mix.Rdata")

pdf("mix_probests_main.pdf",width=7,height=3)
par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(prob1.est.geom2,pwl.prob1.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(prob2.est.geom2,pwl.prob2.est.geom2.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(prob3.est.geom3,pwl.prob3.est.geom3.ss4,
        names=c("Par","PWL"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
dev.off()

pdf("mix_probests_supp.pdf",width=9,height=6)
par(mfrow=c(3,1),pty="m",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
boxplot(pwl.prob1.est.geom2.ss1,pwl.prob1.est.geom2.ss2,pwl.prob1.est.geom2.ss3,pwl.prob1.est.geom2.ss4,
        pwl.prob1.est.geom2.ss5,pwl.prob1.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[8,10]×[0.01,3]")
abline(h=TP1,col="red")
boxplot(pwl.prob2.est.geom2.ss1,pwl.prob2.est.geom2.ss2,pwl.prob2.est.geom2.ss3,pwl.prob2.est.geom2.ss4,
        pwl.prob2.est.geom2.ss5,pwl.prob2.est.geom2.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[5,7]×[0.01,3]")
abline(h=TP2,col="red")
boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
        pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
        names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
        main="[8,10]×[2,4]×[0.01,3]", outline=T)
abline(h=TP3,col="red")
dev.off()

#############################################################################
#############################################################################

rm(list=ls())
load("simstudy1/SimStudy3d_alog1.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3

load("simstudy2/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3

load("simstudy3/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss3 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3

load("simstudy4/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3

load("simstudy5/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3

load("simstudy6/pwlinexp_SimStudy3d_alog1.Rdata")
# pars.ss4 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_alog1.Rdata")

alog1.rmse.par.B1 = sqrt(mean((log(prob1.est.geom2[prob1.est.geom2>0])-log(TP1))^2,na.rm=T))
alog1.rmse.par.B2 = sqrt(mean((log(prob2.est.geom2[prob2.est.geom2>0])-log(TP2))^2,na.rm=T))
alog1.rmse.par.B3 = sqrt(mean((log(prob3.est.geom3[prob3.est.geom3>0])-log(TP3))^2,na.rm=T))

alog1.rmse.pwl1.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss1[pwl.prob1.est.geom2.ss1>0])-log(TP1))^2,na.rm=T))
alog1.rmse.pwl1.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss1[pwl.prob2.est.geom2.ss1>0])-log(TP2))^2,na.rm=T))
alog1.rmse.pwl1.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss1[pwl.prob3.est.geom3.ss1>0])-log(TP3))^2,na.rm=T))

alog1.rmse.pwl2.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss2[pwl.prob1.est.geom2.ss2>0])-log(TP1))^2,na.rm=T))
alog1.rmse.pwl2.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss2[pwl.prob2.est.geom2.ss2>0])-log(TP2))^2,na.rm=T))
alog1.rmse.pwl2.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss2[pwl.prob3.est.geom3.ss2>0])-log(TP3))^2,na.rm=T))

alog1.rmse.pwl3.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss3[pwl.prob1.est.geom2.ss3>0])-log(TP1))^2,na.rm=T))
alog1.rmse.pwl3.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss3[pwl.prob2.est.geom2.ss3>0])-log(TP2))^2,na.rm=T))
alog1.rmse.pwl3.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss3[pwl.prob3.est.geom3.ss3>0])-log(TP3))^2,na.rm=T))

alog1.rmse.pwl4.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss4[pwl.prob1.est.geom2.ss4>0])-log(TP1))^2,na.rm=T))
alog1.rmse.pwl4.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss4[pwl.prob2.est.geom2.ss4>0])-log(TP2))^2,na.rm=T))
alog1.rmse.pwl4.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss4[pwl.prob3.est.geom3.ss4>0])-log(TP3))^2,na.rm=T))

alog1.rmse.pwl5.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss5[pwl.prob1.est.geom2.ss5>0])-log(TP1))^2, na.rm=T))
alog1.rmse.pwl5.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss5[pwl.prob2.est.geom2.ss5>0])-log(TP2))^2, na.rm=T))
alog1.rmse.pwl5.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss5[pwl.prob3.est.geom3.ss5>0])-log(TP3))^2, na.rm=T))

alog1.rmse.pwl6.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss6[pwl.prob1.est.geom2.ss6>0])-log(TP1))^2, na.rm=T))
alog1.rmse.pwl6.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss6[pwl.prob2.est.geom2.ss6>0])-log(TP2))^2, na.rm=T))
alog1.rmse.pwl6.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss6[pwl.prob3.est.geom3.ss6>0])-log(TP3))^2, na.rm=T))

load("simstudy1/SimStudy3d_alog2.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3

load("simstudy2/pwlinexp_SimStudy3d_alog2.Rdata")
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3

load("simstudy3/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3

load("simstudy4/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3

load("simstudy5/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3

load("simstudy6/pwlinexp_SimStudy3d_alog2.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_alog2.Rdata")


alog2.rmse.par.B1 = sqrt(mean((log(prob1.est.geom2[prob1.est.geom2>0])-log(TP1))^2,na.rm=T))
alog2.rmse.par.B2 = sqrt(mean((log(prob2.est.geom2[prob2.est.geom2>0])-log(TP2))^2,na.rm=T))
alog2.rmse.par.B3 = sqrt(mean((log(prob3.est.geom3[prob3.est.geom3>0])-log(TP3))^2,na.rm=T))

alog2.rmse.pwl1.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss1[pwl.prob1.est.geom2.ss1>0])-log(TP1))^2,na.rm=T))
alog2.rmse.pwl1.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss1[pwl.prob2.est.geom2.ss1>0])-log(TP2))^2,na.rm=T))
alog2.rmse.pwl1.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss1[pwl.prob3.est.geom3.ss1>0])-log(TP3))^2,na.rm=T))

alog2.rmse.pwl2.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss2[pwl.prob1.est.geom2.ss2>0])-log(TP1))^2,na.rm=T))
alog2.rmse.pwl2.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss2[pwl.prob2.est.geom2.ss2>0])-log(TP2))^2,na.rm=T))
alog2.rmse.pwl2.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss2[pwl.prob3.est.geom3.ss2>0])-log(TP3))^2,na.rm=T))

alog2.rmse.pwl3.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss3[pwl.prob1.est.geom2.ss3>0])-log(TP1))^2,na.rm=T))
alog2.rmse.pwl3.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss3[pwl.prob2.est.geom2.ss3>0])-log(TP2))^2,na.rm=T))
alog2.rmse.pwl3.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss3[pwl.prob3.est.geom3.ss3>0])-log(TP3))^2,na.rm=T))

alog2.rmse.pwl4.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss4[pwl.prob1.est.geom2.ss4>0])-log(TP1))^2,na.rm=T))
alog2.rmse.pwl4.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss4[pwl.prob2.est.geom2.ss4>0])-log(TP2))^2,na.rm=T))
alog2.rmse.pwl4.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss4[pwl.prob3.est.geom3.ss4>0])-log(TP3))^2,na.rm=T))

alog2.rmse.pwl5.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss5[pwl.prob1.est.geom2.ss5>0])-log(TP1))^2, na.rm=T))
alog2.rmse.pwl5.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss5[pwl.prob2.est.geom2.ss5>0])-log(TP2))^2, na.rm=T))
alog2.rmse.pwl5.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss5[pwl.prob3.est.geom3.ss5>0])-log(TP3))^2, na.rm=T))

alog2.rmse.pwl6.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss6[pwl.prob1.est.geom2.ss6>0])-log(TP1))^2, na.rm=T))
alog2.rmse.pwl6.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss6[pwl.prob2.est.geom2.ss6>0])-log(TP2))^2, na.rm=T))
alog2.rmse.pwl6.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss6[pwl.prob3.est.geom3.ss6>0])-log(TP3))^2, na.rm=T))

load("simstudy1/SimStudy3d_mix.Rdata")
# pars.ss1 = pars
pwl.prob1.est.geom.ss1   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss1  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss1  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss1   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss1  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss1  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss1   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss1  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss1  = pwl.prob3.est.geom3

load("simstudy2/pwlinexp_SimStudy3d_mix.Rdata")
pwl.prob1.est.geom.ss2   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss2  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss2  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss2   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss2  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss2  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss2   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss2  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss2  = pwl.prob3.est.geom3

load("simstudy3/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss3   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss3  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss3  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss3   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss3  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss3  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss3   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss3  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss3  = pwl.prob3.est.geom3

load("simstudy4/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss4   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss4  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss4  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss4   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss4  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss4  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss4   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss4  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss4  = pwl.prob3.est.geom3

load("simstudy5/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss5   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss5  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss5  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss5   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss5  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss5  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss5   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss5  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss5  = pwl.prob3.est.geom3

load("simstudy6/pwlinexp_SimStudy3d_mix.Rdata")
# pars.ss2 = pars
pwl.prob1.est.geom.ss6   = pwl.prob1.est.geom
pwl.prob1.est.geom2.ss6  = pwl.prob1.est.geom2
pwl.prob1.est.geom3.ss6  = pwl.prob1.est.geom3
pwl.prob2.est.geom.ss6   = pwl.prob2.est.geom
pwl.prob2.est.geom2.ss6  = pwl.prob2.est.geom2
pwl.prob2.est.geom3.ss6  = pwl.prob2.est.geom3
pwl.prob3.est.geom.ss6   = pwl.prob3.est.geom
pwl.prob3.est.geom2.ss6  = pwl.prob3.est.geom2
pwl.prob3.est.geom3.ss6  = pwl.prob3.est.geom3

# Probability estimation boxplots
load("../../../pw_lin_gauge/SimStudy_3d/par/simstudypar2/SimStudy3d_mix.Rdata")

mix.rmse.par.B1 = sqrt(mean((log(prob1.est.geom2[prob1.est.geom2>0])-log(TP1))^2,na.rm=T))
mix.rmse.par.B2 = sqrt(mean((log(prob2.est.geom2[prob2.est.geom2>0])-log(TP2))^2,na.rm=T))
mix.rmse.par.B3 = sqrt(mean((log(prob3.est.geom3[prob3.est.geom3>0])-log(TP3))^2,na.rm=T))

mix.rmse.pwl1.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss1[pwl.prob1.est.geom2.ss1>0])-log(TP1))^2,na.rm=T))
mix.rmse.pwl1.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss1[pwl.prob2.est.geom2.ss1>0])-log(TP2))^2,na.rm=T))
mix.rmse.pwl1.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss1[pwl.prob3.est.geom3.ss1>0])-log(TP3))^2,na.rm=T))

mix.rmse.pwl2.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss2[pwl.prob1.est.geom2.ss2>0])-log(TP1))^2,na.rm=T))
mix.rmse.pwl2.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss2[pwl.prob2.est.geom2.ss2>0])-log(TP2))^2,na.rm=T))
mix.rmse.pwl2.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss2[pwl.prob3.est.geom3.ss2>0])-log(TP3))^2,na.rm=T))

mix.rmse.pwl3.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss3[pwl.prob1.est.geom2.ss3>0])-log(TP1))^2,na.rm=T))
mix.rmse.pwl3.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss3[pwl.prob2.est.geom2.ss3>0])-log(TP2))^2,na.rm=T))
mix.rmse.pwl3.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss3[pwl.prob3.est.geom3.ss3>0])-log(TP3))^2,na.rm=T))

mix.rmse.pwl4.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss4[pwl.prob1.est.geom2.ss4>0])-log(TP1))^2,na.rm=T))
mix.rmse.pwl4.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss4[pwl.prob2.est.geom2.ss4>0])-log(TP2))^2,na.rm=T))
mix.rmse.pwl4.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss4[pwl.prob3.est.geom3.ss4>0])-log(TP3))^2,na.rm=T))

mix.rmse.pwl5.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss5[pwl.prob1.est.geom2.ss5>0])-log(TP1))^2, na.rm=T))
mix.rmse.pwl5.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss5[pwl.prob2.est.geom2.ss5>0])-log(TP2))^2, na.rm=T))
mix.rmse.pwl5.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss5[pwl.prob3.est.geom3.ss5>0])-log(TP3))^2, na.rm=T))

mix.rmse.pwl6.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss6[pwl.prob1.est.geom2.ss6>0])-log(TP1))^2, na.rm=T))
mix.rmse.pwl6.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss6[pwl.prob2.est.geom2.ss6>0])-log(TP2))^2, na.rm=T))
mix.rmse.pwl6.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss6[pwl.prob3.est.geom3.ss6>0])-log(TP3))^2, na.rm=T))


all.rmse = matrix(c(alog1.rmse.par.B1,alog2.rmse.par.B1,mix.rmse.par.B1,
                    alog1.rmse.pwl1.B1,alog2.rmse.pwl1.B1,mix.rmse.pwl1.B1,
                    alog1.rmse.pwl2.B1,alog2.rmse.pwl2.B1,mix.rmse.pwl2.B1,
                    alog1.rmse.pwl3.B1,alog2.rmse.pwl3.B1,mix.rmse.pwl3.B1,
                    alog1.rmse.pwl4.B1,alog2.rmse.pwl4.B1,mix.rmse.pwl4.B1,
                    alog1.rmse.pwl5.B1,alog2.rmse.pwl5.B1,mix.rmse.pwl5.B1,
                    alog1.rmse.pwl6.B1,alog2.rmse.pwl6.B1,mix.rmse.pwl6.B1,
                    alog1.rmse.par.B2,alog2.rmse.par.B2,mix.rmse.par.B2,
                    alog1.rmse.pwl1.B2,alog2.rmse.pwl1.B2,mix.rmse.pwl1.B2,
                    alog1.rmse.pwl2.B2,alog2.rmse.pwl2.B2,mix.rmse.pwl2.B2,
                    alog1.rmse.pwl3.B2,alog2.rmse.pwl3.B2,mix.rmse.pwl3.B2,
                    alog1.rmse.pwl4.B2,alog2.rmse.pwl4.B2,mix.rmse.pwl4.B2,
                    alog1.rmse.pwl5.B2,alog2.rmse.pwl5.B2,mix.rmse.pwl5.B2,
                    alog1.rmse.pwl6.B2,alog2.rmse.pwl6.B2,mix.rmse.pwl6.B2,
                    alog1.rmse.par.B3,alog2.rmse.par.B3,mix.rmse.par.B3,
                    alog1.rmse.pwl1.B3,alog2.rmse.pwl1.B3,mix.rmse.pwl1.B3,
                    alog1.rmse.pwl2.B3,alog2.rmse.pwl2.B3,mix.rmse.pwl2.B3,
                    alog1.rmse.pwl3.B3,alog2.rmse.pwl3.B3,mix.rmse.pwl3.B3,
                    alog1.rmse.pwl4.B3,alog2.rmse.pwl4.B3,mix.rmse.pwl4.B3,
                    alog1.rmse.pwl5.B3,alog2.rmse.pwl5.B3,mix.rmse.pwl5.B3,
                    alog1.rmse.pwl6.B3,alog2.rmse.pwl6.B3,mix.rmse.pwl6.B3),
                  ncol=3,byrow=T)
all.rmse
library(xtable)
xtable(all.rmse,digits=4)
