# d=4 Delauney triangulation on air pollution
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

library(PWLExtremes)

###########################

load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_4d_urban.Rdata")

idx.groups = lapply(c(1:4), function(x) combn(4,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

# load("d4pollution_quants_tau070.RData")
# qr = qr.res$qr
# thresh.r = as.numeric(qr$r.tau.wpts)
# w.mesh = qr$wpts
# u.vals.lst = lapply(idx.groups,function(idx){
#   u.min = pexp(max(sapply(1:length(thresh.r),function(w.idx){
#     w.vec = w.mesh[w.idx,]
#     r0w.val = thresh.r[w.idx]
#     return(r0w.val*min(w.vec[idx],na.rm=T))
#   }),na.rm=T))
#   return(seq(u.min,1,length.out=100))
# })
load("u_vals_d4pollution.RData")

gfun = function(w,par,locs=par.locs){
  gfun.pwl(x=w,par=par,ref.angles=locs)
}

load("parlocs2.RData")
par.locs=par.locs.2
# load("R_fit_13_refit.RData")
# init=R.fit.13.refit$init.val

bs.res.pwl = list(list(NULL))
num.bs.reps = 500
for(rep in c(1:num.bs.reps)){
  set.seed(1000+rep)
  
  message(paste("bootstrap rep:",rep))
  
  # block bootstrap
  # acf(ds.exp.4d)
  block.size = 7
  n.idx.sample = nrow(ds.exp.4d)
  idx.init = sample(c(1:nrow(ds.exp.4d)),ceil(nrow(ds.exp.4d)/block.size),replace=T)
  idx = lapply(idx.init,function(ii) ii+c(0:(block.size-1)))
  idx = unlist(idx)
  idx = idx[idx<=nrow(ds.exp.4d)]
  # idx = sample(c(1:nrow(ds.exp.4d)),nrow(ds.exp.4d),replace=T)  # regular bootstrap
  ds.bs = ds.exp.4d[idx,]
  # acf(ds.bs)
  r.bs<-apply(ds.bs,1,sum)
  w.bs<-ds.bs/r.bs
  
  tau=0.70
  qr.bs=PWLExtremes::radial.quants.L1.KDE(r=r.bs,w=w.bs,tau=tau,mesh.eval=FALSE,bww=0.075,ker="Gaussian")#QR.4d.L1.KDE(r=r.bs,w=w.bs,tau=tau,mesh.eval=F)
  excind.bs<-r.bs>qr.bs$r0w
  rexc.bs<-r.bs[excind.bs]
  wexc.bs<-w.bs[excind.bs,]
  r0wexc.bs<-qr.bs$r0w[excind.bs]
  na.ind.bs<-which(is.na(rexc.bs))
  # na.ind
  if(length(na.ind.bs)>0){
    rexc.bs<-rexc.bs[-na.ind.bs]
    wexc.bs<-wexc.bs[-na.ind.bs,]
    r0wexc.bs<-r0wexc.bs[-na.ind.bs]
  }
  
  init = KDE.quant.eval(wpts=par.locs,r=r.bs,w=w.bs,tau=tau,bww=0.075,ker="Gaussian")
  init = init/max(par.locs[,1] * init, na.rm = TRUE)
  
  # parametric and PWL fit on boostrap dataset
  skip_to_next <- FALSE
  tryCatch({
    R.fit.24.bs = fit.pwlin(r=rexc.bs,r0w=r0wexc.bs,w=wexc.bs,locs=par.locs,init.val=init,
                  pen.const=1, 
                  method="BFGS",bound.fit=TRUE)
    W.fit.bs = fit.pwlin(r=rexc.bs,r0w=r0wexc.bs,w=wexc.bs,locs=par.locs,init.val=init/init[1],
                      pen.const=20,method="BFGS",fW.fit=TRUE)
    xstar.pwl.bs = sim.joint(nsim=50000,k.vals=1,shape=4,
                             par.locs=par.locs,par.locs.W=par.locs,
                             par=R.fit.24.bs$mle, fW.par=W.fit.bs$fW.mle,
                             bww=0.075, tau=tau, r=r.bs, w=w.bs)[[1]]
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    message("Error in fitting to bootsrap data.")
    bs.res.pwl[[rep]] = list(NA)
    next
  } else {
    gauge.est.vals.w.bs = gfun.pwl(w.bs,par=R.fit.24.bs$mle,ref.angles=par.locs)
    T.vals = seq(10,1000,by=10)
    est.prop.exc.bs = sapply(T.vals,function(T.val){
      return.radii = qgamma(1-(1/T.val), shape=4, rate=gauge.est.vals.w.bs)
      return(mean(r.bs>return.radii))
    })
    est.log.T.bs = log(1/est.prop.exc.bs)
    
    emp.res = list(list(NULL))
    pwl.res = list(list(NULL))
    
    col.idx=1
    for(which.cols in idx.groups){
      message(which.cols)
      u.vals = u.vals.lst[[col.idx]]
      emp.res[[col.idx]] = sapply(u.vals, function(u) mean(apply(ds.bs[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u))  # assuming the margins are approximately standard exponential
      pwl.res[[col.idx]] = mean(excind.bs,na.rm=T)*sapply(u.vals, function(u) mean(apply(xstar.pwl.bs[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
      col.idx = col.idx + 1
    }
    
    bs.res.pwl[[rep]] = list(pwl.bs.pars=list(
      R.fit.24.bs=R.fit.24.bs$mle
    ),
    pwl.res=list(
      emp.res=emp.res,
      pwl.res=pwl.res
    ),est.log.T.bs = est.log.T.bs)
    save(bs.res.pwl,file=paste0("d4pollution_chi_BSfits.Rdata"))
  }
}
