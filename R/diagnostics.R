# qq.3d<-function(r,w.adj.angles,r0w,shape=3,par,par.locs,quantilefn=qexp,plot=TRUE,...){
#   
#   # QQ plot for exceedance radii
#   
#   g.vals = sapply(w.adj.angles, function(lst){
#     which.angles = lst$idx.locs
#     if(any(is.na(lst$w)) | all(which.angles==0)){  # if angle is not in the positive orthant OR if angle lies on line
#       return(NA)
#     } else {
#       return(gfun.simple.d3pwlin(xyz=lst$w,par=par[which.angles],par.locs=par.locs[which.angles,]))
#     }
#   })
#   
#   num<-pgamma(r,shape=shape,rate=g.vals,lower.tail=F)
#   den<-pgamma(r0w,shape=shape,rate=g.vals,lower.tail=F)
#   
#   nr<-length(r)
#   emp.vals = quantilefn(c(1:nr)/(nr+1))
#   est.vals = quantilefn(sort(1-num/den))
#   
#   if(plot){
#     plot(emp.vals,
#          est.vals,
#          #xlab="empirical",ylab="model",
#          xlab=NA,ylab=NA,
#          pch=20,cex=0.8,...)
#     abline(a=0,b=1)
#     
#     Ulow<-sapply(1:nr,function(i){qbeta(0.025,i,nr+1-i)})
#     Uup<-sapply(1:nr,function(i){qbeta(0.975,i,nr+1-i)})
#     
#     lines(quantilefn(c(1:nr)/(nr+1)),quantilefn(Ulow),lty=2)
#     lines(quantilefn(c(1:nr)/(nr+1)),quantilefn(Uup),lty=2)
#   } 
#   
#   # return a score, MAE of distance
#   return(mean(abs(emp.vals-est.vals)))
# }

# PP and QQ scores

pp.score = function(r, w, r0w, par, gfun, par.locs, n.dims){
  rate = gfun(w, par = par, ref.angles=par.locs)
  num <- pgamma(r, shape = n.dims, rate = rate, lower.tail = F)
  den <- pgamma(r0w, shape = n.dims, rate = rate, lower.tail = F)
  nr <- length(r)
  score = mean(abs(c(1:nr)/(nr + 1) - sort(1 - num/den)))
  return(score)
}

qq.score = function(r, w, r0w, par, gfun, par.locs, n.dims){
  rate = gfun(w, par = par, ref.angles=par.locs)
  num <- pgamma(r, shape = n.dims, rate = rate, lower.tail = F)
  den <- pgamma(r0w, shape = n.dims, rate = rate, lower.tail = F)
  nr <- length(r)
  score = mean(abs(qexp(c(1:nr)/(nr + 1)) - qexp(sort(1 - num/den))))
  return(score)
}

# return level scores

return.level.score = function(r,w,gfun,par,par.locs,tau,n.dims){
  gauge.est.vals.w = gfun(w,par=par,ref.angles=par.locs)
  T.vals = seq(1/(1-tau),1000,length.out=100)
  est.prop.exc = sapply(T.vals,function(T.val){
    return.radii = qgamma(1-(1/T.val), shape=n.dims, rate=gauge.est.vals.w)
    return(mean(r>return.radii))
  })
  true.prop.exc = 1/T.vals
  
  score = mean(abs(est.prop.exc - true.prop.exc))
  return(score)
}


########################

## PP and QQ plots for truncated Gamma in d=5, not available in the geometricMVE package
ppdiag.5d = function (r, w, r0w, par, gfun, add.gauge = FALSE, gauge1, gauge2, 
          par1ind, par2ind) {
  if (dim(w)[2] == 4) {
    w <- cbind(w, 1 - w[, 1] - w[, 2] - w[, 3] - w[, 4])
  }
  if (!add.gauge) {
    num <- pgamma(r, shape = par[1], rate = apply(w, 1, gfun, par = par[2:length(par)]), lower.tail = F)
    den <- pgamma(r0w, shape = par[1], rate = apply(w, 1, gfun, par = par[2:length(par)]), lower.tail = F)
  }
  else {
    wgrid <- create.wgrid.4d(20)
    ms <- additivegauge.scaling.4d(gauge1 = gauge1, gauge2 = gauge2, 
                                   weight = par[length(par)], par1 = par[par1ind], 
                                   par2 = par[par2ind], wgrid = wgrid)
    num <- pgamma(r, shape = par[1], rate = apply(w, 1, 
                                                  additivegauge.rescale.4d, par1 = par[par1ind], par2 = par[par2ind], 
                                                  ms = ms, gauge1 = gauge1, gauge2 = gauge2, weight = par[length(par)]), 
                  lower.tail = F)
    den <- pgamma(r0w, shape = par[1], rate = apply(w, 1, 
                                                    additivegauge.rescale.4d, par1 = par[par1ind], par2 = par[par2ind], 
                                                    ms = ms, gauge1 = gauge1, gauge2 = gauge2, weight = par[length(par)]), 
                  lower.tail = F)
  }
  nr <- length(r)
  plot(c(1:nr)/(nr + 1), sort(1 - num/den), xlab = "empirical", 
       ylab = "model", pch = 20, cex = 0.8)
  abline(a = 0, b = 1)
  Ulow <- sapply(1:nr, function(i) {
    qbeta(0.025, i, nr + 1 - i)
  })
  Uup <- sapply(1:nr, function(i) {
    qbeta(0.975, i, nr + 1 - i)
  })
  lines(c(1:nr)/(nr + 1), Ulow, lty = 2)
  lines(c(1:nr)/(nr + 1), Uup, lty = 2)
}

qqdiag.5d = function (r, w, r0w, par, gfun, add.gauge = FALSE, gauge1, gauge2, 
                      par1ind, par2ind, quantilefn = qexp) 
{
  if (dim(w)[2] == 4) {
    w <- cbind(w, 1 - w[, 1] - w[, 2] - w[, 3] - w[, 4])
  }
  if (!add.gauge) {
    num <- pgamma(r, shape = par[1], rate = apply(w, 1, gfun, par = par[2:length(par)]), lower.tail = F)
    den <- pgamma(r0w, shape = par[1], rate = apply(w, 1, gfun, par = par[2:length(par)]), lower.tail = F)
  }
  else {
    wgrid <- create.wgrid.4d(20)
    ms <- additivegauge.scaling.4d(gauge1 = gauge1, gauge2 = gauge2, 
                                   weight = par[length(par)], par1 = par[par1ind], 
                                   par2 = par[par2ind], wgrid = wgrid)
    num <- pgamma(r, shape = par[1], rate = apply(w, 1, 
                                                  additivegauge.rescale.4d, par1 = par[par1ind], par2 = par[par2ind], 
                                                  ms = ms, gauge1 = gauge1, gauge2 = gauge2, weight = par[length(par)]), 
                  lower.tail = F)
    den <- pgamma(r0w, shape = par[1], rate = apply(w, 1, 
                                                    additivegauge.rescale.4d, par1 = par[par1ind], par2 = par[par2ind], 
                                                    ms = ms, gauge1 = gauge1, gauge2 = gauge2, weight = par[length(par)]), 
                  lower.tail = F)
  }
  nr <- length(r)
  plot(quantilefn(c(1:nr)/(nr + 1)), quantilefn(sort(1 - num/den)), 
       xlab = "empirical", ylab = "model", pch = 20, cex = 0.8)
  abline(a = 0, b = 1)
  Ulow <- sapply(1:nr, function(i) {
    qbeta(0.025, i, nr + 1 - i)
  })
  Uup <- sapply(1:nr, function(i) {
    qbeta(0.975, i, nr + 1 - i)
  })
  lines(quantilefn(c(1:nr)/(nr + 1)), quantilefn(Ulow), lty = 2)
  lines(quantilefn(c(1:nr)/(nr + 1)), quantilefn(Uup), lty = 2)
}
