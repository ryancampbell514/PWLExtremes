# d=3 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(lattice)
library(scales)
library(geometricMVE)

setwd("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/d3")

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))
library(PWLExtremes)

# the perspective  at which we want to view rgl figures
usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

###############################################################################

# OBTAIN CHI VALUES

load("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/ds_3d_urban.Rdata")
r<-apply(ds.exp.3d,1,sum)
w<-ds.exp.3d/r

# chi plots - empirical vs. geometricMVE vs. pwlin
idx.groups = lapply(c(1:3), function(x) combn(3,x))
idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
idx.groups = unlist(idx.groups, recursive = F)
idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]

# get u values of length 100
# u.min = pexp((1/3)*KDE.quant.eval.3d(wpts = rep(1/3,3), r=r, w=w))
# u.vals = seq(u.min,1,length.out=100)
load("u_vals_d3pollution.RData")

bs.res.emp = list(list(NULL))
num.bs.reps = 500
for(rep in c(1:num.bs.reps)){
  set.seed(1000+rep)
  
  message(paste("bootstrap rep:",rep))
  
  # block bootstrap
  # acf(ds.exp.3d)
  block.size = 7
  n.idx.sample = nrow(ds.exp.3d)
  idx.init = sample(c(1:nrow(ds.exp.3d)),ceil(nrow(ds.exp.3d)/block.size),replace=T)
  idx = lapply(idx.init,function(ii) ii+c(0:(block.size-1)))
  idx = unlist(idx)
  idx = idx[idx<=nrow(ds.exp.3d)]
  # idx = sample(c(1:nrow(ds.exp.3d)),nrow(ds.exp.3d),replace=T)   # regular bootstrap
  ds.bs = ds.exp.3d[idx,]
  # acf(ds.bs)
  
  emp.res = list(list(NULL))
  
  col.idx=1
  for(which.cols in idx.groups){
    message(which.cols)
    
    u.vals = u.vals.lst[[col.idx]]
    
    # empirical
    emp.res[[col.idx]] = sapply(u.vals, function(u) {
      mean(apply(ds.bs[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u)
    })
    col.idx = col.idx + 1
  }
  
  bs.res.emp[[rep]] = emp.res
  save(bs.res.emp,file=paste0("d3pollution_chi_bs_emp.Rdata"))

}
