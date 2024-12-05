# 2-dimensional code

# importance weights:
iweights.2d.pwl = function (k, r0w, w, gfun, shape, par){
  if (k < 1) {
    warning("k is below 1, adjusted to k = 1.")
    k = 1
  }
  rate <- sapply(w, gfun, par = par)
  return(pgamma(k * r0w, shape = shape, rate = rate, lower.tail = F)/
           pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
}

# simulate from the empirical angular distirbution:
sim.2d.cond = function (w, r0w, k = 1, nsim, shape=2, par, gfun, marg="pos") {
  if (k != 1) {
    iw <- iweights.2d.pwl(k = k, r0w = r0w, w = w, gfun = gfun, shape = shape, par = par)
    star.ind <- sample(1:length(w), size = nsim, replace = T, prob = iw)
    wstar <- w[star.ind]
    r0w_star <- c(k * r0w)[star.ind]
  }
  else {
    star.ind <- sample(1:length(w), size = nsim, replace = T)
    wstar <- w[star.ind]
    r0w_star <- r0w[star.ind]
  }
  rate0 <- sapply(w, gfun, par = par)
  rate <- rate0[star.ind]
  rstar <- qgamma(1 - runif(nsim) * pgamma(r0w_star, shape = shape, 
                                           rate = rate, lower.tail = F), 
                  shape = shape, rate = rate)
  
  if(marg=="pos"){
    xstar <- cbind(rstar * wstar, rstar * (1 - wstar))
  } else if(marg=="Rd"){
    xstar <- pol2cart.L1Rd(r=rstar,w=wstar)
  }
  return(xstar)
}

# simulate from the proposed angular density
fW.mcmc.g.2d<-function(niter,nburn,bpar1=1,bpar2=1,thin,g,return.acc.rate=FALSE,...){
  ## bpar1,bpar2 -> beta parameters from proposal distribution, (default is uniform)
  ## ... -> arguments to pass onto the gauge function
  
  # w<-runif(1)
  w = rbeta(1,bpar1,bpar2)
  draws<-numeric(niter)
  it<-1
  acc=0
  while(it<=niter){
    wcan<-rbeta(1,bpar1,bpar2)
    accn<-dbeta(w,bpar1,bpar2)*(g(wcan,...)^(-2))
    accd<-dbeta(wcan,bpar1,bpar2)*(g(w,...)^(-2))
    if(runif(1)<accn/accd){w<-wcan; acc = acc+1}
    draws[it]<-w
    it<-it+1
  }
  print(paste("W MCMC acc. rate:",round(acc/niter,4)))
  if(return.acc.rate){
    return(list(acc.rate=round(acc/niter,4),
                sample=draws[-(1:nburn)][seq(1,(niter-nburn),by=thin)]))
  } else {
    return(draws[-(1:nburn)][seq(1,(niter-nburn),by=thin)])
  }
}

sim.2d.joint = function(nsim,k.vals=1,gfun,shape=2,par,fW.par,par.locs,r,w,
                        tau=0.95,bww=0.05,bwr=0.05,marg="pos"){
  ## k.vals  -> a single, or a vector of k values
  ## par     -> parameters of R|W
  ## fW.par  -> parameters of W
  
  # fit the proposal Beta distribution
  if(marg=="Rd"){
    w.beta = (w+2)/4
    beta.nll = function(pars){-sum(dbeta(x=w.beta,shape1=pars[1],shape2=pars[2],log=T))}
  } else {
    beta.nll = function(pars){-sum(dbeta(x=w,shape1=pars[1],shape2=pars[2],log=T))}
  }
  beta.mle = optim(par=rep(1,2),fn=beta.nll)$par
  
  nthin=2
  wstar = fW.mcmc.g.2d(niter=(nthin*nsim)+1000,nburn=1000,
                       bpar1=beta.mle[1],bpar2=beta.mle[2],thin=nthin,
                       g=gfun,par=fW.par,locs=par.locs)
  
  if(marg=="Rd"){
    wstar = wstar*4 -2
  }
  
  # now, evaluate r0w at the sampled angles
  r.tau.wstar = sapply(wstar,function(wstar.i){
    weightsw<-dnorm(w,mean=as.numeric(wstar.i),sd=bww)
    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,30))
    return(ur$root)
  })
  
  sims = lapply(k.vals, function(k){
    if(k==1){
      rate <- sapply(wstar, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(k*r.tau.wstar, shape = shape, 
                                               rate = rate, lower.tail = F), 
                      shape = shape, rate = rate)
      # return(rstar*cbind(wstar,1-wstar))
      if(marg=="pos"){
        xstar <- cbind(rstar * wstar, rstar * (1 - wstar))
      } else if(marg=="Rd"){
        xstar <- pol2cart.L1Rd(r=rstar,w=wstar)
      }
      return(xstar)
    } else if (k>1){
      iw <- iweights.2d.pwl(k=k,r0w=r.tau.wstar,w=wstar,gfun=gfun,shape=shape,par=par)
      star.ind <- sample(1:length(wstar),size=nsim,replace=T,prob=iw)
      wstar2 <- wstar[star.ind]
      r.tau.wstar2 <- c(k*r.tau.wstar)[star.ind]
      
      rate <- sapply(wstar2, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(r.tau.wstar2, shape = shape, 
                                               rate = rate, lower.tail = F), 
                      shape = shape, rate = rate)
      # return(rstar*cbind(wstar2,1-wstar2))
      if(marg=="pos"){
        xstar <- cbind(rstar * wstar2, rstar * (1 - wstar2))
      } else if(marg=="Rd"){
        xstar <- pol2cart.L1Rd(r=rstar,w=wstar2)
      }
      return(xstar)
    }
  })
  return(sims)
}


##############################################################################
##############################################################################
# general d-dimensions

iweights.pwl = function (k, r0w, w, gfun, shape, par){
  # WARNING: assumed fixed shape parameter of d
  
  if (k < 1) {
    warning("k is below 1, adjusted to k = 1.")
    k = 1
  }
  rate=apply(w, 1, gfun, par = par)
  return(pgamma(k * r0w, shape = shape, rate = rate, lower.tail = F)/
           pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
}

sim.cond = function (w, r0w, k = 1, nsim, shape=NULL, par, gfun) {
  
  if(is.null(shape)){
    shape=dim(w)[2]
  }
  
  # WARNING: assumed fixed shape parameter
  if (k != 1) {
    iw <- iweights.pwl(k = k, r0w = r0w, w = w, gfun = gfun, par = par)
    star.ind <- sample(1:nrow(w), size = nsim, replace = T, prob = iw)
    wstar <- w[star.ind,]
    r0w_star <- c(k * r0w)[star.ind]
  }
  else {
    star.ind <- sample(1:nrow(w), size = nsim, replace = T)
    wstar <- w[star.ind,]
    r0w_star <- r0w[star.ind]
  }
  # rate0 <- pwlin.g.vals.3d(which.adj.angles.3d(angles=w,locs=locs),par,locs)
  rate0=apply(w, 1, gfun, par = par)
  rate <- rate0[star.ind]
  rstar <- qgamma(1 - runif(nsim) * pgamma(r0w_star, shape = shape, 
                                           rate = rate, lower.tail = F), 
                  shape = shape, rate = rate)
  xstar = rstar * wstar
  return(xstar)
}

# sample using angle information
fW.mcmc.g<-function(W.dim,niter,nburn,alpha,thin,g,return.acc.rate=FALSE,...){
  ## alpha -> vector, parameters of the Dirichlet distribution (default -> uniform on simplex)
  ## ...   -> arguments to pass onto the gauge function
  
  if(is.null(alpha)){
    alpha=rep(1,W.dim)
  }
  
  # initialize
  w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
  draws<-matrix(ncol=W.dim,nrow=niter)
  it<-1
  acc=0
  while(it<=niter){
    wcan<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
    
    accn<-LaplacesDemon::ddirichlet(w,alpha)*(g(wcan,...)^(-W.dim))
    accd<-LaplacesDemon::ddirichlet(wcan,alpha)*(g(w,...)^(-W.dim))
    
    if(runif(1)<accn/accd){w<-wcan; acc=acc+1}
    
    draws[it,]<-w
    it<-it+1
  }
  print(paste("W MCMC acc. rate:",round(acc/niter,4)))
  if(return.acc.rate){
    return(list(acc.rate=round(acc/niter,4),
                sample=draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),]))
  } else {
    return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
  }
}

sim.joint = function(nsim,k.vals=1,shape=NULL,par,fW.par,par.locs,par.locs.W=NULL,
                     r,w,tau=0.95,bww=0.05,bwr=0.05){
  
  if(is.null(par.locs.W)){
    par.locs.W=par.locs
  }
  
  # get proposal alphas
  is.error <- FALSE
  tryCatch({
    dirichlet.nll = function(pars){-sum(LaplacesDemon::ddirichlet(x=w,alpha=pars,log=T))}
  },error=function(e){
    is.error <<- TRUE 
  })
  
  W.dim=dim(w)[2]
  if(is.null(shape)){
    shape=W.dim
  }
  
  if(is.error) {
    alpha=rep(1,W.dim)  # uniform over the simplex
  } else {
    alpha = optim(par=rep(1,W.dim),fn=dirichlet.nll)$par
  }
  # dirichlet.nll = function(pars){-sum(LaplacesDemon::ddirichlet(x=w,alpha=pars,log=T))}
  # alpha = optim(par=rep(1,3),fn=dirichlet.nll)$par
  # alpha=c(1,1,1)  # uniform over the simplex
  
  # simulate exceedance angles
  nthin=4
  wstar = fW.mcmc.g(W.dim=W.dim,niter=(nthin*nsim)+1000,nburn=1000,thin=nthin,
                       alpha=alpha,
                       g=function(w,par=fW.par,locs=par.locs.W){gfun.pwl(x=w,par=par,ref.angles=locs)},
                    par=fW.par,locs=par.locs.W)
  
  # now, evaluate r0w at the sampled angles
  r.tau.wstar = apply(wstar,1,function(wstar.i){
    weightsw<-dmvnorm(w[,-W.dim],mean=as.numeric(wstar.i[-W.dim]),sigma=bww^2*diag(W.dim-1))
    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,30))
    return(ur$root)
  })
  
  gfun = function(w,par=par,locs=par.locs){gfun.pwl(x=w,par=par,ref.angles=locs)}
  
  sims = lapply(k.vals, function(k){
    if(k==1){
      rate <- apply(wstar, 1, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(k*r.tau.wstar, shape = shape, 
                                               rate = rate, lower.tail = F),
                      shape = shape, rate = rate)
      return(rstar*wstar)
    } else if (k>1){
      iw <- iweights.pwl(k=k,r0w=r.tau.wstar,w=wstar,gfun=gfun,shape=shape,par=par)
      star.ind <- sample(1:nrow(wstar),size=nsim,replace=T,prob=iw)
      wstar2 <- wstar[star.ind,]
      r.tau.wstar2 <- c(k*r.tau.wstar)[star.ind]
      
      rate <- apply(wstar2, 1, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(r.tau.wstar2, shape = shape, 
                                               rate = rate, lower.tail = F), 
                      shape = shape, rate = rate)
      return(rstar*wstar2)
    }
  })
  return(sims)
}

