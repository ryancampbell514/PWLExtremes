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

###############################################################################

# OBTAIN CHI VALUES FOR MODEL FITS

# load("d3_pen_model.RData")
# par.locs = d3.pollution.fits$par.locs
# mle = d3.pollution.fits$joint.model.fit.pen.bounded$mle

load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")
# load("d3pollution_fits_tau95.RData")
# load("d3pollution_fits_finermesh.RData")
# r<-apply(ds.exp.3d,1,sum)
# w<-ds.exp.3d/r
# tau=0.9
# qr=d3.pollution.fits$qr
# r0w = qr$r0w
# excind<-r>r0w
# rexc<-r[excind]
# wexc<-w[excind,]
# r0w<-r0w[excind]

idx.groups = lapply(c(1:3), function(x) combn(3,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

# get u values of length 100
load("u_vals_d3pollution.RData")

gfun = function(w,par,locs=par.locs){
  adj.angles = which.adj.angles.3d(angles=w, locs)
  pwlin.g.vals.3d(adj.angles,par,locs)
}

###############################################################################

# OBTAIN THE BOOTSTRAP FITS AND CHI VALUES
# block.size = 15  # from acf of original data
# num.blocks = ceiling(dim(ds.exp.3d)[1]/block.size)

par.locs = seq(0,1,by=1/6)   # newly-suggested
par.locs = as.matrix(expand.grid(par.locs,par.locs))
par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
par.locs = par.locs[apply(par.locs,1,function(w) !any(w<0)),]
par.locs[,3] = ifelse(par.locs[,3]<0.001,0,par.locs[,3])
nnodes = nrow(par.locs)

# bs.res.pwl = list(list(NULL))
load("d3pollution_chi_BSfits.Rdata")
num.bs.reps = 500
start.idx = 80
for(rep in c(start.idx:num.bs.reps)){
  set.seed(1000+rep)
  
  message(paste("bootstrap rep:",rep))
  
  # block bootstrap
  block.size = 7
  n.idx.sample = nrow(ds.exp.3d)
  idx.init = sample(c(1:nrow(ds.exp.3d)),ceil(nrow(ds.exp.3d)/block.size),replace=T)
  idx = lapply(idx.init,function(ii) ii+c(0:(block.size-1)))
  idx = unlist(idx)
  idx = idx[idx<=nrow(ds.exp.3d)]
  # idx = sample(c(1:nrow(ds.exp.3d)),nrow(ds.exp.3d),replace=T)   # regular bootstrap
  ds.bs = ds.exp.3d[idx,]
  r.bs<-apply(ds.bs,1,sum)
  w.bs<-ds.bs/r.bs
  
  tau=0.95
  qr.bs=PWLExtremes::radial.quants.L1.KDE(r=r.bs,w=w.bs,tau=tau,mesh.eval=FALSE,bww=0.075,ker="Gaussian")
  excind.bs<-r.bs>qr.bs$r0w
  rexc.bs<-r.bs[excind.bs]
  wexc.bs<-w.bs[excind.bs,]
  r0wexc.bs<-qr.bs$r0w[excind.bs]
  na.ind.bs<-which(is.na(rexc.bs))
  # na.ind
  if(length(na.ind.bs)>0){
    rexc.bs<-rexc.bs[-na.ind.bs]
    wexc.bs<-wexc.bs[-na.ind.bs,]
    r0wexc.bs<-r0wexc.bs[-na.ind.bs]}
  
  init = KDE.quant.eval(wpts=par.locs,r=r.bs,w=w.bs,tau=tau,bww=0.075,ker="Gaussian")
  init = init/max(par.locs[,1] * init, na.rm = TRUE)
  
  # parametric and PWL fit on boostrap dataset
  skip_to_next <- FALSE
  tryCatch({
    R.fit.24.bs = fit.pwlin(r=rexc.bs,r0w=r0wexc.bs,w=wexc.bs,locs=par.locs,pen.const=1,method="BFGS",init.val=init,bound.fit=TRUE)#bound.fit.3d(r=rexc.bs,w=wexc.bs,r0w=r0w.bs,par.locs,pen.const=best.pen.const,method="BFGS")
    W.fit.34.bs = fit.pwlin(r=rexc.bs,r0w=r0wexc.bs,w=wexc.bs,locs=par.locs,pen.const=20,method="BFGS",joint.fit=F,fW.fit=T,init.val=init/init[1])
    
    xstar.pwl.bs.4 = sim.joint(nsim=50000,k.vals=1,shape=3,
              par.locs=par.locs,par.locs.W=par.locs,
              par=R.fit.24.bs$mle, fW.par=W.fit.34.bs$fW.mle,
              bww=0.075, tau=tau, r=r.bs, w=w.bs)[[1]]
  },error=function(e){
    skip_to_next <<- TRUE
  })
  if(skip_to_next) {
    message("Error in fitting to bootsrap data.")
    bs.res.pwl[[rep]] = list(NA)
    next
  } else {
    # return-level plot estimation
    gauge.est.vals.w.bs = gfun.pwl(w.bs,par=R.fit.24.bs$mle,ref.angles=par.locs)
    T.vals = seq(21,1000,by=10)
    est.prop.exc.bs = sapply(T.vals,function(T.val){
      return.radii = qgamma(1-(1/T.val), shape=3, rate=gauge.est.vals.w.bs)
      return(mean(r.bs>return.radii))
    })
    est.log.T.bs = log(1/est.prop.exc.bs)
    
    # chi estimation
    pwl.res = emp.res = list(list(NULL))
    
    col.idx=1
    for(which.cols in idx.groups){
      message(which.cols)
      u.vals = u.vals.lst[[col.idx]]
      emp.res[[col.idx]] = sapply(u.vals, function(u) {
        mean(apply(ds.bs[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u)
      })
      pwl.res[[col.idx]] = mean(excind.bs,na.rm=T)*sapply(u.vals, function(u) mean(apply(xstar.pwl.bs.4[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
      col.idx = col.idx + 1
    }
    
    bs.res.pwl[[rep]] = list(pwl.bs.pars=list(W.fit.34.bs=W.fit.34.bs$fW.mle,
                                              R.fit.24.bs=R.fit.24.bs$mle),
                         bs.idx=idx,
                         excind.bs=excind.bs,
                         pwl.res=pwl.res,
                         emp.res=emp.res,
                         est.log.T.bs = est.log.T.bs)
    save(bs.res.pwl,file=paste0("d3pollution_chi_BSfits.Rdata"))
  }
}


