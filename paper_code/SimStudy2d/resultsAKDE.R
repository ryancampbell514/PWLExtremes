# Compile and plot d=2 simstudy results to show best model
rm(list=ls())
# library(geometricMVE)
library(scales)
library(this.path)

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# library(PWLExtremes)

setwd("/mnt/usb-TOSHIBA_External_USB_3.0_20200304013644F-0:0-part1/lancaster/pw_lin_gauge/SimStudy_2d")

fn.dir = "../../geometricMVE/R"
invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

for(cop in c("log_strongdep","log_weakdep","gauss","invlog")){
  
  load(paste0("pwl_gradpen_newpen_fixedlambda/simstudy1/pwlinexp_SimStudy2d_",cop,".Rdata"))
  old.pwl.prob1.est.geom.ss1  = pwl.prob1.est.geom
  old.pwl.prob1.est.geom2.ss1 = pwl.prob1.est.geom2
  old.pwl.prob2.est.geom.ss1  = pwl.prob2.est.geom
  old.pwl.prob2.est.geom2.ss1 = pwl.prob2.est.geom2
  old.pwl.prob3.est.geom.ss1  = pwl.prob3.est.geom
  old.pwl.prob3.est.geom2.ss1 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    old.pwl.prob1.est.geom3.ss1 = pwl.prob1.est.geom3
  } else {
    old.pwl.prob3.est.geom3.ss1 = pwl.prob3.est.geom3
  }
  
  load(paste0("pwl_gitcode_AKDEthresh/simstudy1/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss1 = pars
  new.pwl.prob1.est.geom.ss1  = pwl.prob1.est.geom
  new.pwl.prob1.est.geom2.ss1 = pwl.prob1.est.geom2
  new.pwl.prob2.est.geom.ss1  = pwl.prob2.est.geom
  new.pwl.prob2.est.geom2.ss1 = pwl.prob2.est.geom2
  new.pwl.prob3.est.geom.ss1  = pwl.prob3.est.geom
  new.pwl.prob3.est.geom2.ss1 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    new.pwl.prob1.est.geom3.ss1 = pwl.prob1.est.geom3
  } else {
    new.pwl.prob3.est.geom3.ss1 = pwl.prob3.est.geom3
  }
  
  load(paste0("pwl_gradpen_newpen_fixedlambda/simstudy2/pwlinexp_SimStudy2d_",cop,".Rdata"))
  old.pwl.prob1.est.geom.ss2  = pwl.prob1.est.geom
  old.pwl.prob1.est.geom2.ss2 = pwl.prob1.est.geom2
  old.pwl.prob2.est.geom.ss2  = pwl.prob2.est.geom
  old.pwl.prob2.est.geom2.ss2 = pwl.prob2.est.geom2
  old.pwl.prob3.est.geom.ss2  = pwl.prob3.est.geom
  old.pwl.prob3.est.geom2.ss2 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    old.pwl.prob1.est.geom3.ss2 = pwl.prob1.est.geom3
  } else {
    old.pwl.prob3.est.geom3.ss2 = pwl.prob3.est.geom3
  }
  
  load(paste0("pwl_gitcode_AKDEthresh/simstudy2/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss2 = pars
  new.pwl.prob1.est.geom.ss2  = pwl.prob1.est.geom
  new.pwl.prob1.est.geom2.ss2 = pwl.prob1.est.geom2
  new.pwl.prob2.est.geom.ss2  = pwl.prob2.est.geom
  new.pwl.prob2.est.geom2.ss2 = pwl.prob2.est.geom2
  new.pwl.prob3.est.geom.ss2  = pwl.prob3.est.geom
  new.pwl.prob3.est.geom2.ss2 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    new.pwl.prob1.est.geom3.ss2 = pwl.prob1.est.geom3
  } else {
    new.pwl.prob3.est.geom3.ss2 = pwl.prob3.est.geom3
  }
  
  load(paste0("/mnt/usb-TOSHIBA_External_USB_3.0_20200304013644F-0:0-part1/lancaster/pw_lin_gauge/SimStudy_2d/par/simstudypar2/SimStudy2d_",cop,".Rdata"))
  
  # library(scales)
  # par(pty="s")
  # plot(x/log(nrow(x)),col=alpha("grey",0.3),pch=20)
  # for(ii in 1:nrow(pars)){
  #   lines(cbind(qr$wpts,1-qr$wpts)/
  #           apply(cbind(qr$wpts,1-qr$wpts),1,gauge_rvad,par=pars[ii,2]),
  #         col=alpha("blue",0.1))
  # }
  
  pdf(paste0("pwl_gitcode_AKDEthresh/SimStudy2d_",cop,"_SS1probests_main_compareAKDE.pdf"),width=7,height=3)
  par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
  
  if(cop=="invlog"){
    boxplot(prob1.est.geom3,old.pwl.prob1.est.geom3.ss1,new.pwl.prob1.est.geom3.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[10,12]",outline=T)
    abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,old.pwl.prob2.est.geom2.ss1,new.pwl.prob2.est.geom2.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom2,old.pwl.prob3.est.geom2.ss1,new.pwl.prob3.est.geom2.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[2,4]")
    abline(h=prob3.true,col="red")
  } else {
    boxplot(prob1.est.geom2,old.pwl.prob1.est.geom2.ss1,new.pwl.prob1.est.geom2.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[10,12]")
    abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,old.pwl.prob2.est.geom2.ss1,new.pwl.prob2.est.geom2.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom3,old.pwl.prob3.est.geom3.ss1,new.pwl.prob3.est.geom3.ss1,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[2,4]",outline=T)
    abline(h=prob3.true,col="red")
  }
  dev.off()
  
  pdf(paste0("pwl_gitcode_AKDEthresh/SimStudy2d_",cop,"_SS2probests_main_compareAKDE.pdf"),width=7,height=3)
  par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)
  
  if(cop=="invlog"){
    boxplot(prob1.est.geom3,old.pwl.prob1.est.geom3.ss2,new.pwl.prob1.est.geom3.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[10,12]",outline=T)
    abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,old.pwl.prob2.est.geom2.ss2,new.pwl.prob2.est.geom2.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom2,old.pwl.prob3.est.geom2.ss2,new.pwl.prob3.est.geom2.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[2,4]")
    abline(h=prob3.true,col="red")
  } else {
    boxplot(prob1.est.geom2,old.pwl.prob1.est.geom2.ss2,new.pwl.prob1.est.geom2.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[10,12]")
    abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,old.pwl.prob2.est.geom2.ss2,new.pwl.prob2.est.geom2.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom3,old.pwl.prob3.est.geom3.ss2,new.pwl.prob3.est.geom3.ss2,
            names=c("Par","PWL","AKDE"),
            main="[10,12]×[2,4]",outline=T)
    abline(h=prob3.true,col="red")
  }
  dev.off()
  
  # plot the gauge functions
  pars.ss1.med = apply(pars.ss1,2,median,na.rm=T)
  pars.ss2.med = apply(pars.ss2,2,median,na.rm=T)

  gfun = function(w,par,locs=par.locs){
    adj.angles = which.adj.angles.2d(angles=w, locs)
    pwlin.g.vals.2d(adj.angles,par,locs)
  }

  gauge.med.vals.ss1 = sapply(ww,gfun,par=pars.ss1.med,locs=seq(0,1,length.out=length(pars.ss1.med)))
  gauge.med.vals.ss2 = sapply(ww,gfun,par=pars.ss2.med,locs=seq(0,1,length.out=length(pars.ss2.med)))

  if(cop %in% c("log_strongdep","log_weakdep")){
    true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_rvad,par=aAD)
  } else if(cop=="gauss"){
    true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_gaussian,par=rho)
  } else if(cop=="invlog"){
    true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_invlogistic,par=aAI)
  } else {
    stop("select a valid copula")
  }

  est.gauge.vals.ss1 = apply(pars.ss1[!apply(pars.ss1,1,function(par) all(is.na(par))),],1,function(gauge.pars){
    sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  })
  est.gauge.vals.ss1 = t(est.gauge.vals.ss1)
  est.gauge.vals.ss2 = apply(pars.ss2[!apply(pars.ss2,1,function(par) all(is.na(par))),],1,function(gauge.pars){
    sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  })
  est.gauge.vals.ss2 = t(est.gauge.vals.ss2)

  pdf(paste0("pwl_gitcode_AKDEthresh/AKDE_SimStudy2d_",cop,"_SS1_unitg.pdf"),width=4,height=4)
  par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  for(row in 1:nrow(est.gauge.vals.ss1)){
    if(row==1){
      plot(cbind(ww,1-ww)/est.gauge.vals.ss1[row,],
           type="l",
           col=scales::alpha("grey",0.3),
           xlim=c(0,1.1),ylim=c(0,1.1),
           xlab=NA,ylab=NA,
           main=NA)
    } else {
      lines(cbind(ww,1-ww)/est.gauge.vals.ss1[row,],
            col=scales::alpha("grey",0.3))
    }
  }
  lines(x=c(0,1), y=c(1,1), lty=3)
  lines(x=c(1,1), y=c(0,1), lty=3)
  lines(cbind(ww,1-ww)/gauge.med.vals.ss1,lwd=2,col="blue")
  lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  dev.off()

  pdf(paste0("pwl_gitcode_AKDEthresh/AKDE_SimStudy2d_",cop,"_SS2_unitg.pdf"),width=4,height=4)
  par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  for(row in 1:nrow(est.gauge.vals.ss2)){
    if(row==1){
      plot(cbind(ww,1-ww)/est.gauge.vals.ss2[row,],
           type="l",
           col=scales::alpha("grey",0.3),
           xlim=c(0,1.1),ylim=c(0,1.1),
           xlab=NA,ylab=NA,
           main=NA)
    } else {
      lines(cbind(ww,1-ww)/est.gauge.vals.ss2[row,],
            col=scales::alpha("grey",0.3))
    }
  }
  lines(x=c(0,1), y=c(1,1), lty=3)
  lines(x=c(1,1), y=c(0,1), lty=3)
  lines(cbind(ww,1-ww)/gauge.med.vals.ss2,lwd=2,col="blue")
  lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  dev.off()
  
}
