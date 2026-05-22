# Compile and plot d=2 simstudy results to show best model
rm(list=ls())
# library(geometricMVE)
library(scales)
library(this.path)

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
# library(PWLExtremes)

setwd(path.join(this.path::here(),"pwl_gitcode_AKDEthresh_simplex"))

##############################################################################

all.rmse = NULL

for(cop in c("log_strongdep","log_weakdep","gauss","invlog")){
  print(cop)
  
  load(paste0("simstudy1/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss1 = pars
  pwl.prob1.est.geom.ss1  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss1 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss1  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss1 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss1  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss1 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss1 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss1 = pwl.prob3.est.geom3
  }
  # plot(x,pch=20,col="grey")
  # lines(qr$r.tau.wpts*cbind(qr$wpts,1-qr$wpts))
  
  load(paste0("simstudy2/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss2 = pars
  pwl.prob1.est.geom.ss2  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss2 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss2  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss2 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss2  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss2 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss2 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss2 = pwl.prob3.est.geom3
  }  
  # plot(x,pch=20,col="grey")
  # lines(qr$r.tau.wpts*cbind(qr$wpts,1-qr$wpts))
  
  load(paste0("simstudy3/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss3 = pars
  fW.pars.ss3 = mle.fW
  wexc.ss3 = wexc
  wexc.samp.ss3 = xstar.sims[[1]] / apply(xstar.sims[[1]],1,sum)
  pwl.prob1.est.geom.ss3  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss3 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss3  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss3 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss3  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss3 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss3 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss3 = pwl.prob3.est.geom3
  }

  load(paste0("simstudy4/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss4 = pars
  pwl.prob1.est.geom.ss4  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss4 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss4  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss4 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss4  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss4 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss4 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss4 = pwl.prob3.est.geom3
  }
  
  load(paste0("simstudy5/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss5 = pars
  fW.pars.ss5 = mle.fW
  wexc.ss5 = wexc
  wexc.samp.ss5 = xstar.sims[[1]] / apply(xstar.sims[[1]],1,sum)
  pwl.prob1.est.geom.ss5  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss5 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss5  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss5 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss5  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss5 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss5 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss5 = pwl.prob3.est.geom3
  }
  
  load(paste0("simstudy6/pwlinexp_SimStudy2d_",cop,".Rdata"))
  pars.ss6 = pars
  fW.pars.ss6 = mle.fW
  wexc.ss6 = wexc
  wexc.samp.ss6 = xstar.sims[[1]] / apply(xstar.sims[[1]],1,sum)
  pwl.prob1.est.geom.ss6  = pwl.prob1.est.geom
  pwl.prob1.est.geom2.ss6 = pwl.prob1.est.geom2
  pwl.prob2.est.geom.ss6  = pwl.prob2.est.geom
  pwl.prob2.est.geom2.ss6 = pwl.prob2.est.geom2
  pwl.prob3.est.geom.ss6  = pwl.prob3.est.geom
  pwl.prob3.est.geom2.ss6 = pwl.prob3.est.geom2
  if(cop=="invlog"){
    pwl.prob1.est.geom3.ss6 = pwl.prob1.est.geom3
  } else {
    pwl.prob3.est.geom3.ss6 = pwl.prob3.est.geom3
  }

  # ##############################################################################
  # 
  # # density of angles
  # nb = ifelse(cop=="gauss",20,30)
  # pdf(paste0("SimStudy2d_",cop,"_fW.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.ss3,breaks=nb,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # # evaluate the density of angles at the median parameters
  # lines(ww,
  #       0.5*(G.vol.2d(fW.pars.ss3,par.locs)^(-1))*sapply(ww,gfun,par=fW.pars.ss3)^(-2),
  #       lwd=2,col="blue")
  # lines(ww,
  #       0.5*(G.vol.2d(fW.pars.ss5,par.locs)^(-1))*sapply(ww,gfun,par=fW.pars.ss5)^(-2),
  #       lwd=2,col="darkgreen",lty=2)
  # lines(ww,
  #       0.5*(G.vol.2d(fW.pars.ss6,par.locs)^(-1))*sapply(ww,gfun,par=fW.pars.ss6)^(-2),
  #       lwd=2,col="red",lty=4)
  # legend("topright",legend=c("SS3/4","SS5","SS6"),lty=c(1,2,4),lwd=c(2,2,2),col=c("blue","darkgreen","red"))
  # dev.off()
  # 
  # pdf(paste0("SimStudy2d_",cop,"_SS3_fW.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.ss3,breaks=nb,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # # evaluate the density of angles at the median parameters
  # lines(ww,
  #       0.5*(G.vol.2d(fW.pars.ss3,par.locs)^(-1))*sapply(ww,gfun,par=fW.pars.ss3)^(-2),
  #       lwd=2,col="blue")
  # dev.off()
  # 
  # # Wexc histograms
  # pdf(paste0("SimStudy2d_",cop,"_emp_Wexc_samp.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.ss3,breaks=nb,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # dev.off()
  # pdf(paste0("SimStudy2d_",cop,"_SS3_Wexc_samp.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.samp.ss3,breaks=30,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # dev.off()
  # pdf(paste0("SimStudy2d_",cop,"_SS5_Wexc_samp.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.samp.ss5,breaks=30,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # dev.off()
  # pdf(paste0("SimStudy2d_",cop,"_SS6_Wexc_samp.pdf"),width=5,height=5)
  # par(mfrow=c(1,1),pty="s")
  # hist(wexc.samp.ss6,breaks=30,freq=F,main=NA,ylim=c(0,6),xlab="w")
  # dev.off()
  # 
  # # plot the gauge functions
  # pars.ss1.med = apply(pars.ss1,2,median,na.rm=T)
  # pars.ss2.med = apply(pars.ss2,2,median,na.rm=T)
  # pars.ss5.med = apply(pars.ss5,2,median,na.rm=T)
  # pars.ss6.med = apply(pars.ss6,2,median,na.rm=T)
  # 
  # gfun = function(w,par,locs=par.locs){
  #   adj.angles = which.adj.angles.2d(angles=w, locs)
  #   pwlin.g.vals.2d(adj.angles,par,locs)
  # }
  # 
  # gauge.med.vals.ss1 = sapply(ww,gfun,par=pars.ss1.med,locs=seq(0,1,length.out=length(pars.ss1.med)))
  # gauge.med.vals.ss2 = sapply(ww,gfun,par=pars.ss2.med,locs=seq(0,1,length.out=length(pars.ss2.med)))
  # gauge.med.vals.ss5 = sapply(ww,gfun,par=pars.ss5.med,locs=seq(0,1,length.out=length(pars.ss5.med)))
  # gauge.med.vals.ss6 = sapply(ww,gfun,par=pars.ss6.med,locs=seq(0,1,length.out=length(pars.ss6.med)))
  # 
  # if(cop %in% c("log_strongdep","log_weakdep")){
  #   true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_rvad,par=aAD)
  # } else if(cop=="gauss"){
  #   true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_gaussian,par=rho)
  # } else if(cop=="invlog"){
  #   true.gauge.vals = apply(cbind(ww,1-ww),1,gauge_invlogistic,par=aAI)
  # } else {
  #   stop("select a valid copula")
  # }
  # 
  # est.gauge.vals.ss1 = apply(pars.ss1[!apply(pars.ss1,1,function(par) all(is.na(par))),],1,function(gauge.pars){
  #   sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  # })
  # est.gauge.vals.ss1 = t(est.gauge.vals.ss1)
  # est.gauge.vals.ss2 = apply(pars.ss2[!apply(pars.ss2,1,function(par) all(is.na(par))),],1,function(gauge.pars){
  #   sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  # })
  # est.gauge.vals.ss2 = t(est.gauge.vals.ss2)
  # est.gauge.vals.ss5 = apply(pars.ss5[!apply(pars.ss5,1,function(par) all(is.na(par))),],1,function(gauge.pars){
  #   sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  # })
  # est.gauge.vals.ss5 = t(est.gauge.vals.ss5)
  # est.gauge.vals.ss6 = apply(pars.ss6[!apply(pars.ss6,1,function(par) all(is.na(par))),],1,function(gauge.pars){
  #   sapply(ww,gfun,par=gauge.pars,locs=seq(0,1,length.out=length(gauge.pars)))
  # })
  # est.gauge.vals.ss6 = t(est.gauge.vals.ss6)
  # 
  # pdf(paste0("SimStudy2d_",cop,"_SS1_unitg.pdf"),width=4,height=4)
  # par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  # for(row in 1:nrow(est.gauge.vals.ss1)){
  #   if(row==1){
  #     plot(cbind(ww,1-ww)/est.gauge.vals.ss1[row,],
  #          type="l",
  #          col=scales::alpha("grey",0.3),
  #          xlim=c(0,1.1),ylim=c(0,1.1),
  #          xlab=NA,ylab=NA,
  #          main=NA)
  #   } else {
  #     lines(cbind(ww,1-ww)/est.gauge.vals.ss1[row,],
  #           col=scales::alpha("grey",0.3))
  #   }
  # }
  # lines(x=c(0,1), y=c(1,1), lty=3)
  # lines(x=c(1,1), y=c(0,1), lty=3)
  # lines(cbind(ww,1-ww)/gauge.med.vals.ss1,lwd=2,col="blue")
  # lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  # dev.off()
  # 
  # pdf(paste0("SimStudy2d_",cop,"_SS2_unitg.pdf"),width=4,height=4)
  # par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  # for(row in 1:nrow(est.gauge.vals.ss2)){
  #   if(row==1){
  #     plot(cbind(ww,1-ww)/est.gauge.vals.ss2[row,],
  #          type="l",
  #          col=scales::alpha("grey",0.3),
  #          xlim=c(0,1.1),ylim=c(0,1.1),
  #          xlab=NA,ylab=NA,
  #          main=NA)
  #   } else {
  #     lines(cbind(ww,1-ww)/est.gauge.vals.ss2[row,],
  #           col=scales::alpha("grey",0.3))
  #   }
  # }
  # lines(x=c(0,1), y=c(1,1), lty=3)
  # lines(x=c(1,1), y=c(0,1), lty=3)
  # lines(cbind(ww,1-ww)/gauge.med.vals.ss2,lwd=2,col="blue")
  # lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  # dev.off()
  # 
  # pdf(paste0("SimStudy2d_",cop,"_SS5_unitg.pdf"),width=4,height=4)
  # par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  # for(row in 1:nrow(est.gauge.vals.ss5)){
  #   if(row==1){
  #     plot(cbind(ww,1-ww)/est.gauge.vals.ss5[row,],
  #          type="l",
  #          col=scales::alpha("grey",0.3),
  #          xlim=c(0,1.1),ylim=c(0,1.1),
  #          xlab=NA,ylab=NA,
  #          main=NA)
  #   } else {
  #     lines(cbind(ww,1-ww)/est.gauge.vals.ss5[row,],
  #           col=scales::alpha("grey",0.3))
  #   }
  # }
  # lines(x=c(0,1), y=c(1,1), lty=3)
  # lines(x=c(1,1), y=c(0,1), lty=3)
  # lines(cbind(ww,1-ww)/gauge.med.vals.ss5,lwd=2,col="blue")
  # lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  # dev.off()
  # 
  # pdf(paste0("SimStudy2d_",cop,"_SS6_unitg.pdf"),width=4,height=4)
  # par(mfrow=c(1,1),pty="s",mar=c(3,1,1,1))
  # for(row in 1:nrow(est.gauge.vals.ss6)){
  #   if(row==1){
  #     plot(cbind(ww,1-ww)/est.gauge.vals.ss6[row,],
  #          type="l",
  #          col=scales::alpha("grey",0.3),
  #          xlim=c(0,1.1),ylim=c(0,1.1),
  #          xlab=NA,ylab=NA,
  #          main=NA)
  #   } else {
  #     lines(cbind(ww,1-ww)/est.gauge.vals.ss6[row,],
  #           col=scales::alpha("grey",0.3))
  #   }
  # }
  # lines(x=c(0,1), y=c(1,1), lty=3)
  # lines(x=c(1,1), y=c(0,1), lty=3)
  # lines(cbind(ww,1-ww)/gauge.med.vals.ss6,lwd=2,col="blue")
  # lines(cbind(ww,1-ww)/true.gauge.vals,lwd=1,lty=2,col="black")
  # dev.off()

  ##############################################################################

  # Probability estimation boxplots
  load(paste0("../../../pw_lin_gauge/SimStudy_2d/par/simstudypar2/SimStudy2d_",cop,".Rdata"))

  # library(scales)
  # par(pty="s")
  # plot(x/log(nrow(x)),col=alpha("grey",0.3),pch=20)
  # for(ii in 1:nrow(pars)){
  #   lines(cbind(qr$wpts,1-qr$wpts)/
  #           apply(cbind(qr$wpts,1-qr$wpts),1,gauge_rvad,par=pars[ii,2]),
  #         col=alpha("blue",0.1))
  # }

  pdf(paste0("SimStudy2d_",cop,"_probests_main.pdf"),width=7,height=3)
  par(mfrow=c(1,3),pty="s",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)

  if(cop=="invlog"){
    boxplot(prob1.est.geom3,pwl.prob1.est.geom3.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[10,12]",outline=T)
    abline(h=prob1.true,col="red")
    # boxplot(prob1.est.geom3,pwl.prob1.est.geom3.ss4,
    #         names=c("Par","PWL"),
    #         main="[10,12]×[10,12] (no outliers)",outline=F)
    # abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,pwl.prob2.est.geom2.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom2,pwl.prob3.est.geom2.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[2,4]")
    abline(h=prob3.true,col="red")
  } else {
    boxplot(prob1.est.geom2,pwl.prob1.est.geom2.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[10,12]")
    abline(h=prob1.true,col="red")
    boxplot(prob2.est.geom2,pwl.prob2.est.geom2.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(prob3.est.geom3,pwl.prob3.est.geom3.ss4,
            names=c("Par","PWL"),
            main="[10,12]×[2,4]",outline=T)
    abline(h=prob3.true,col="red")
    # boxplot(prob3.est.geom3,pwl.prob3.est.geom3.ss4,
    #         names=c("Par","PWL"),
    #         main="[10,12]×[2,4] (no outliers)",outline=F)
    # abline(h=prob3.true,col="red")
  }
  dev.off()

  pdf(paste0("SimStudy2d_",cop,"_probests_supp.pdf"),width=9,height=9)
  par(mfrow=c(4,1),pty="m",mar=c(2.5,2.5,1.5,1),cex.main=1.5,cex.axis=1.5)

  if(cop=="invlog"){
    boxplot(pwl.prob1.est.geom3.ss1,pwl.prob1.est.geom3.ss2,
            pwl.prob1.est.geom3.ss3,pwl.prob1.est.geom3.ss4,
            pwl.prob1.est.geom3.ss5,pwl.prob1.est.geom3.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[10,12]",outline=T)
    abline(h=prob1.true,col="red")
    boxplot(pwl.prob1.est.geom3.ss1,pwl.prob1.est.geom3.ss2,
            pwl.prob1.est.geom3.ss3,pwl.prob1.est.geom3.ss4,
            pwl.prob1.est.geom3.ss5,pwl.prob1.est.geom3.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[10,12] (no outliers)",outline=F)
    abline(h=prob1.true,col="red")
    boxplot(pwl.prob2.est.geom2.ss1,pwl.prob2.est.geom2.ss2,
            pwl.prob2.est.geom2.ss3,pwl.prob2.est.geom2.ss4,
            pwl.prob2.est.geom2.ss5,pwl.prob2.est.geom2.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(pwl.prob3.est.geom2.ss1,pwl.prob3.est.geom2.ss2,
            pwl.prob3.est.geom2.ss3,pwl.prob3.est.geom2.ss4,
            pwl.prob3.est.geom2.ss5,pwl.prob3.est.geom2.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[2,4]")
    abline(h=prob3.true,col="red")
  } else {
    boxplot(pwl.prob1.est.geom2.ss1,pwl.prob1.est.geom2.ss2,
            pwl.prob1.est.geom2.ss3,pwl.prob1.est.geom2.ss4,
            pwl.prob1.est.geom2.ss5,pwl.prob1.est.geom2.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[10,12]")
    abline(h=prob1.true,col="red")
    boxplot(pwl.prob2.est.geom2.ss1,pwl.prob2.est.geom2.ss2,
            pwl.prob2.est.geom2.ss3,pwl.prob2.est.geom2.ss4,
            pwl.prob2.est.geom2.ss5,pwl.prob2.est.geom2.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[6,8]")
    abline(h=prob2.true,col="red")
    boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,
            pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
            pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[2,4]",outline=T)
    abline(h=prob3.true,col="red")
    boxplot(pwl.prob3.est.geom3.ss1,pwl.prob3.est.geom3.ss2,
            pwl.prob3.est.geom3.ss3,pwl.prob3.est.geom3.ss4,
            pwl.prob3.est.geom3.ss5,pwl.prob3.est.geom3.ss6,
            names=c("SS1","SS2","SS3","SS4","SS5","SS6"),
            main="[10,12]×[2,4] (no outliers)",outline=F)
    abline(h=prob3.true,col="red")
  }
  dev.off()
  
  if(cop=="invlog"){
    rmse.par.B1 = sqrt(mean((log(prob1.est.geom3)-log(prob1.true))^2))
    rmse.par.B2 = sqrt(mean((log(prob2.est.geom2)-log(prob2.true))^2))
    rmse.par.B3 = sqrt(mean((log(prob3.est.geom2)-log(prob3.true))^2))
      
    rmse.pwl1.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss1)-log(prob1.true))^2))
    rmse.pwl1.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss1)-log(prob2.true))^2))
    rmse.pwl1.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss1)-log(prob3.true))^2))
    
    rmse.pwl2.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss2)-log(prob1.true))^2))
    rmse.pwl2.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss2)-log(prob2.true))^2))
    rmse.pwl2.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss2)-log(prob3.true))^2))
    
    rmse.pwl3.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss3)-log(prob1.true))^2))
    rmse.pwl3.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss3)-log(prob2.true))^2))
    rmse.pwl3.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss3)-log(prob3.true))^2))
    
    rmse.pwl4.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss4)-log(prob1.true))^2))
    rmse.pwl4.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss4)-log(prob2.true))^2))
    rmse.pwl4.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss4)-log(prob3.true))^2))
    
    rmse.pwl5.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss5)-log(prob1.true))^2))
    rmse.pwl5.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss5)-log(prob2.true))^2))
    rmse.pwl5.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss5)-log(prob3.true))^2))
    
    rmse.pwl6.B1 = sqrt(mean((log(pwl.prob1.est.geom3.ss6)-log(prob1.true))^2))
    rmse.pwl6.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss6)-log(prob2.true))^2))
    rmse.pwl6.B3 = sqrt(mean((log(pwl.prob3.est.geom2.ss6)-log(prob3.true))^2))
      
  } else {
    rmse.par.B1 = sqrt(mean((log(prob1.est.geom2)-log(prob1.true))^2))
    rmse.par.B2 = sqrt(mean((log(prob2.est.geom2)-log(prob2.true))^2))
    rmse.par.B3 = sqrt(mean((log(prob3.est.geom3)-log(prob3.true))^2))
    
    rmse.pwl1.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss1)-log(prob1.true))^2))
    rmse.pwl1.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss1)-log(prob2.true))^2))
    rmse.pwl1.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss1[pwl.prob3.est.geom3.ss1>0])-log(prob3.true))^2))
    
    rmse.pwl2.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss2)-log(prob1.true))^2))
    rmse.pwl2.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss2)-log(prob2.true))^2))
    rmse.pwl2.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss2[pwl.prob3.est.geom3.ss2>0])-log(prob3.true))^2))
    
    rmse.pwl3.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss3)-log(prob1.true))^2))
    rmse.pwl3.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss3)-log(prob2.true))^2))
    rmse.pwl3.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss3[pwl.prob3.est.geom3.ss3>0])-log(prob3.true))^2))
    
    rmse.pwl4.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss4)-log(prob1.true))^2))
    rmse.pwl4.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss4)-log(prob2.true))^2))
    rmse.pwl4.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss4[pwl.prob3.est.geom3.ss4>0])-log(prob3.true))^2))
    
    rmse.pwl5.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss5)-log(prob1.true))^2))
    rmse.pwl5.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss5)-log(prob2.true))^2))
    rmse.pwl5.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss5[pwl.prob3.est.geom3.ss5>0])-log(prob3.true))^2))
    
    rmse.pwl6.B1 = sqrt(mean((log(pwl.prob1.est.geom2.ss6)-log(prob1.true))^2))
    rmse.pwl6.B2 = sqrt(mean((log(pwl.prob2.est.geom2.ss6)-log(prob2.true))^2))
    rmse.pwl6.B3 = sqrt(mean((log(pwl.prob3.est.geom3.ss6[pwl.prob3.est.geom3.ss6>0])-log(prob3.true))^2))
    
  }
  
  all.rmse = cbind(all.rmse,
                   c(rmse.par.B1,rmse.pwl1.B1,rmse.pwl2.B1,rmse.pwl3.B1,rmse.pwl4.B1,rmse.pwl5.B1,rmse.pwl6.B1,
                     rmse.par.B2,rmse.pwl1.B2,rmse.pwl2.B2,rmse.pwl3.B2,rmse.pwl4.B2,rmse.pwl5.B2,rmse.pwl6.B2,
                     rmse.par.B3,rmse.pwl1.B3,rmse.pwl2.B3,rmse.pwl3.B3,rmse.pwl4.B3,rmse.pwl5.B3,rmse.pwl6.B3))

  # # Histogram of penalty parameters
  # png(paste0("SimStudy2d_",cop,"_penstrength.png"),width=500,height=500)
  # par(mfrow=c(1,2),pty="m",mar=c(2.5,2.5,1.5,1))
  # hist(pen.consts.ss5,
  #      breaks=length(unique(pen.consts.ss5)),
  #      xlab="pen. strength",
  #      main="simstudy5")
  # hist(pen.consts.ss6,
  #      breaks=length(unique(pen.consts.ss6)),
  #      xlab="pen. strength",
  #      main="simstudy6")
  # dev.off()
}

all.rmse
library(xtable)
xtable(all.rmse,digits=4)

