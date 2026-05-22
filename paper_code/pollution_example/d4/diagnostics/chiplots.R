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


# chi plots

idx.groups = lapply(c(1:4), function(x) combn(4,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

thresh.r = as.numeric(qr$r.tau.wpts)
w.mesh = qr$wpts
# u.vals.lst = lapply(idx.groups,function(idx){
#   u.min = pexp(max(sapply(1:length(thresh.r),function(w.idx){
#     w.vec = w.mesh[w.idx,]
#     r0w.val = thresh.r[w.idx]
#     return(r0w.val*min(w.vec[idx],na.rm=T))
#   }),na.rm=T))
#   return(seq(u.min,1,length.out=100))
# })
# save(u.vals.lst,file="u_vals_d4pollution.RData")
load("u_vals_d4pollution.RData")

# chi.vals.emp = list(list(NULL))
# gfun = function(w,par,locs=par.locs){
#   gfun.pwl(w,par,locs)
# }
xstar.sims = xstar.sims.joint#xstar.sims.joint#sim.cond(w=wexc,r0w=r0w,nsim=50000,k=1,gfun=gfun,par=mle)
col.names = as.character(c(1:4))
chi.vals = list(list(NULL))
grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = u.vals.lst[[grp.no]]
  
  # empirical
  chi.u.emp = sapply(u.vals, function(u) mean(apply(ds.exp.4d[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u))  # assuming the margins are approximately standard exponential
  
  # pwl
  chi.u.est.pwl = mean(excind,na.rm=T)*sapply(u.vals, function(u) mean(apply(xstar.sims[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
  
  col.names = as.character(c(1:4))
  col.names = col.names[which.cols]
  name = paste(col.names, collapse = "")
  
  chi.vals[[grp.no]] = list(u.vals=u.vals,
                            chi.u.emp=chi.u.emp,
                            chi.u.est.pwl=chi.u.est.pwl,
                            name=name)
  
  # ylab = bquote(.(rlang::sym("chi"))[.(name)](u))
  # 
  # # pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_chi_plot_",name,".pdf"),height=5,width=5)
  # par(mfrow=c(1,1),pty="s",mar=c(4,4.5,1,4.5),cex.axis=1.5,cex.lab=1.5)
  # plot(u.vals,chi.u.emp,type="l",lwd=1,
  #      xlab="u",ylab=ylab,
  #      ylim=c(0,1.0))
  # lines(u.vals,chi.u.est.pwl,lty=2,col="blue",lwd=2)
  # # dev.off()
  
  grp.no=grp.no+1
}

grp.no=1
for(which.cols in idx.groups){
  message(which.cols)
  
  u.vals = chi.vals[[grp.no]]$u.vals
  chi.u.emp = chi.vals[[grp.no]]$chi.u.emp
  emp.CI = do.call(rbind,lapply(bs.res.pwl,function(lst){
    lst$pwl.res$emp.res[[grp.no]]
  }))
  emp.CI.low = apply(emp.CI,2,quantile,p=0.025,na.rm=T)
  emp.CI.upp = apply(emp.CI,2,quantile,p=0.975,na.rm=T)
  chi.u.est.pwl = chi.vals[[grp.no]]$chi.u.est.pwl
  pwl.CI =  do.call(rbind,lapply(bs.res.pwl,function(lst){
    lst$pwl.res$pwl.res[[grp.no]]
  }))
  pwl.CI.low = apply(pwl.CI,2,quantile,p=0.025,na.rm=T)
  pwl.CI.upp = apply(pwl.CI,2,quantile,p=0.975,na.rm=T)
  
  col.names = as.character(c(1:4))
  col.names = col.names[which.cols]
  name = chi.vals[[grp.no]]$name
  
  ylab = bquote(.(rlang::sym("chi"))[.(name)](u))
  
  pdf(paste0("~/Dropbox/phd_research/pw_lin_gauge/paper/figs_EVApaper/d4pollution_chi_plot_",name,".pdf"),height=5,width=5)
  par(mfrow=c(1,1),pty="s",mar=c(4,4.5,1,4.5),cex.axis=1.5,cex.lab=1.5)
  plot(u.vals,chi.u.emp,type="l",lwd=1,
       xlab="u",ylab=ylab,
       ylim=c(0,1.0))
  CI.x = c(u.vals,rev(u.vals))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(emp.CI.low,rev(emp.CI.upp))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=alpha("gray",0.3),border=NA)
  lines(u.vals,chi.u.est.pwl,lty=2,col="blue",lwd=2)
  CI.x = c(u.vals,rev(u.vals))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(pwl.CI.low,rev(pwl.CI.upp))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=alpha("blue",0.25),border=NA)
  dev.off()
  
  grp.no=grp.no+1
}
