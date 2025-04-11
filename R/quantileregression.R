# UNIVARIATE KERNELS...

# The Gaussian kernel density and distribution functions
Gaussian.ker.pdf = function(x,mean,sd){dnorm(x,mean=mean,sd=sd)}
mv.Gaussian.ker.pdf = function(x,mean,sd){
  dd = length(mean)
  dmvnorm(x,mean=mean,sigma=(sd^2)*diag(dd))
}
Gaussian.ker.cdf = function(x,mean,sd){pnorm(x,mean=mean,sd=sd)}

# The Epanechnikov kernel density and distribution functions
Epanechnikov.ker.pdf = function(x,mean,sd){
  u = (x-mean)/sd
  val = pmax((1/sd)*0.75*(1-u^2),0)
  return(val)
}
Epanechnikov.ker.cdf = function(x,mean,sd){
  u = (x-mean)/sd
  val = 0.75*((u*(1-((1/3)*u^2))) + (2/3))#(1/sd)*0.75*((u*(1-((1/3)*u^2))) + (2/3))
  cond1 = u < -1
  cond2 = u > 1
  val[cond1] = 0
  val[cond2] = 1
  return(val)
}

# UNIVARIATE COPULA DENSITY (for wweights)
mv.ker.pdf = function(x,mean,sd,ker.pdf){

  x = x[,-ncol(x)]
  mean = mean[-length(mean)]

  val = lapply(1:nrow(x),function(row) prod(ker.pdf(x=x[row,],mean=mean,sd=sd),na.rm=T))
  val = unlist(val)

  return(val)
}

###############################

get_bww.2d = function(r=r,w=w,tau=tau,bwr=bwr,ker.pdf=ker.pdf,ker.cdf=ker.cdf){

  bww.vals = seq(0,0.2,length.out=20)[-1]  # don't think we can have a bandwidth of zero

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  check.scores = rep(NA,length(bww.vals))

  iter=1
  for(bww.idx in c(1:length(bww.vals))){

    bww = bww.vals[bww.idx]
    # pen.const.W = pen.consts.W[pen.const.idx]

    # lik.scores.R13 = rep(NA,k.folds)

    # for(K in 1:k.folds){
    K = 1
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K]
    r.eval = r[grps==K]
    w.eval = w[grps==K]

    is_error <- FALSE
    tryCatch({
      r0w.eval = sapply(w.eval,function(wpts.i){
        weightsw<-ker.pdf(w.fitting,mean=wpts.i,sd=bww)
        pos.weights = weightsw>0
        weightsw = weightsw[pos.weights]
        ccdf<-function(rc){
          ker.vals = ker.cdf(rc,mean=r.fitting,sd=bwr)[pos.weights]
          num = weightsw*ker.vals
          denom = weightsw
          # print(num)
          # print(denom)
          sum(num,na.rm=T)/sum(denom,na.rm=T)
        }
        dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
        ur<-uniroot(dummy,interval = c(0,30))
        return(ur$root)
      })

      check.scores[bww.idx] = mean(checkfn(quant=tau,r.eval-r0w.eval),na.rm=T)

    },error=function(e){
      is_error <<- TRUE
    })
    if(is_error) {
      # print(pen.const)
      check.scores[bww.idx] = 1e10
      next
    }
  }
  return(bww.vals[which.min(check.scores)])
}

get_bww = function(r=r,w=w,tau=tau,bwr=bwr,ker.pdf=ker.pdf,ker.cdf=ker.cdf){

  bww.vals = seq(0,0.2,length.out=20)[-1]

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  check.scores = rep(NA,length(bww.vals))

  num.cols = dim(w)[2]

  iter=1
  for(bww.idx in c(1:length(bww.vals))){

    bww = bww.vals[bww.idx]
    # pen.const.W = pen.consts.W[pen.const.idx]

    # lik.scores.R13 = rep(NA,k.folds)

    # for(K in 1:k.folds){
    K = 1
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K,]
    r.eval = r[grps==K]
    w.eval = w[grps==K,]

    # is_error <- FALSE
    # tryCatch({
      r0w.eval = apply(w.eval,1,function(wpts.i){
        if(sum(wpts.i)>1){
          return(NA)
        } else{
          weightsw = dmvnorm(w.fitting[,-num.cols],mean=wpts.i[-num.cols],sigma=(bww^2)*diag(num.cols-1))

          ccdf<-function(rc){
            mean(weightsw*Gaussian.ker.cdf(rc,mean=r.fitting,sd=bwr))/mean(weightsw)
          }
          dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
          ur<-uniroot(dummy,interval = c(0,30))
          return(ur$root)
        }
      })

      check.scores[bww.idx] = mean(checkfn(quant=tau,r.eval-r0w.eval),na.rm=T)

    # },error=function(e){
    #   is_error <<- TRUE
    # })
    # if(is_error) {
    #   # print(pen.const)
    #   check.scores[bww.idx] = 1e10
    #   next
    # }
  }
  return(bww.vals[which.min(check.scores)])
}


# kernel density estimation of the cdf when angles are defined by the L1-norm, then invert it
radial.quants.L1.KDE.2d = function(r,w,tau=0.95,bww=0.05,bwr=0.05,
                                   ker.pdf=Gaussian.ker.pdf,ker.cdf=Gaussian.ker.cdf,
                                   mesh.eval=TRUE,n.mesh=200){
  # r, w              -> vectors
  # bww, bwr          -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)
  # n.mesh            -> mesh for wpts
  # ker.pdf, ker.cdf  -> kernel pdf and cdf functions

  if(is.null(bww)){
    message("searching for optimal angular bandwidth for KDE threshold...")
    bww = get_bww.2d(r=r,w=w,tau=tau,bwr=bwr,ker.pdf=ker.pdf,ker.cdf=ker.cdf)
    message(paste("Fitting threshold at angular bandwidth",bww))
  }

  r0w = sapply(w, function(ww){
    weightsw<-ker.pdf(w,mean=ww,sd=bww)

    ccdf<-function(rc){
      mean(weightsw*ker.cdf(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,30))
    return(ur$root)
  })

  if(mesh.eval){
    wpts<-seq(0,1,len=n.mesh)
    r.tau.wpts = sapply(wpts,function(wpts.i){
      weightsw<-ker.pdf(w,mean=wpts.i,sd=bww)
      pos.weights = weightsw>0
      weightsw = weightsw[pos.weights]
      ccdf<-function(rc){
        ker.vals = ker.cdf(rc,mean=r,sd=bwr)[pos.weights]
        num = weightsw*ker.vals
        denom = weightsw
        # print(num)
        # print(denom)
        sum(num,na.rm=T)/sum(denom,na.rm=T)
      }
      dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
      ur<-uniroot(dummy,interval = c(0,30))
      return(ur$root)
    })

    return(list(wpts=wpts,
                r.tau.wpts=r.tau.wpts,
                r0w=r0w))
  } else {
    return(list(r0w=r0w))
  }
}



## General d-dimensional version of QR
radial.quants.L1.KDE = function(r,w,tau=0.95,bww=0.05,bwr=0.05,n.mesh=30,
                        mesh.eval=TRUE, ker.pdf=mv.Gaussian.ker.pdf, ker.cdf=Gaussian.ker.cdf){
  require(mvtnorm)
  # r, w      -> vector and matrix, angles defined by the L1 norm
  # bww, bwr  -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)

  if(is.null(bww)){
    message("searching for optimal angular bandwidth for KDE threshold...")
    bww = get_bww(r=r,w=w,tau=tau,bwr=bwr,ker.pdf=ker.pdf,ker.cdf=ker.cdf)
    message(paste("Fitting threshold at angular bandwidth",bww))
  }

  num.cols = dim(w)[2]

  r0w = apply(w, 1, function(ww){

    # weightsw = dmvnorm(w[,-num.cols],mean=ww[-num.cols],sigma=(bww^2)*diag(num.cols-1))
    weightsw = ker.pdf(x=w[,-num.cols], mean=ww[-num.cols], sd=bww)

    ccdf<-function(rc){
      mean(weightsw*Gaussian.ker.cdf(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,30))
    return(ur$root)
  })

  if(mesh.eval){
    wpts<-expand.grid(replicate(num.cols-1, seq(0,1,len=n.mesh), simplify = FALSE))#expand.grid(seq(0,1,len=n.mesh),seq(0,1,len=n.mesh),seq(0,1,len=n.mesh),seq(0,1,len=n.mesh))
    r.tau.wpts = apply(wpts,1,function(wpts.i){
      if(sum(wpts.i)>1){
        return(NA)
      } else{
        weightsw = dmvnorm(w[,-num.cols],mean=wpts.i[-num.cols],sigma=(bww^2)*diag(num.cols-1))

        ccdf<-function(rc){
          mean(weightsw*Gaussian.ker.cdf(rc,mean=r,sd=bwr))/mean(weightsw)
        }
        dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
        ur<-uniroot(dummy,interval = c(0,30))
        return(ur$root)
      }
    })

    return(list(wpts=cbind(wpts,1-apply(wpts,1,sum)),
                r.tau.wpts=r.tau.wpts,
                r0w=r0w))
  } else {
    return(list(r0w=r0w))
  }
}

#########################

# the (univariate) check function:
checkfn = function(quant,inp){
  inp*(quant - ifelse(inp<0,1,0))
}

quant.score.2d.checkfn = function(smoothness.lvl,w,r,tau=0.95,k.folds=5,method=c("empirical","KDE"),...){
  # Perform K-fold cross-validation on the bww parameter
  # r,w -> vectors, defined by L1 norm
  # tau -> the quantile
  # r0w -> the tau thrshold value at w

  if(min(w)<0 | max(w)>1){
    stop("angles should be on the L1 Simplex.")
  }

  # w.mesh = seq(0,1,length.out=nbins)

  grps = sample(1:k.folds,length(r),replace=T)

  score = mean(sapply(1:k.folds,function(K){
    # split data
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K]
    r.eval = r[grps==K]
    w.eval = w[grps==K]

    # fit the quantile estimates
    if(method=="KDE"){
      r0w.eval = KDE.quant.eval.2d(wpts=w.eval,r=r.fitting,w=w.fitting,tau=tau,bww=smoothness.lvl,bwr=0.05,...)
    } else if (method=="empirical"){
      r0w.eval = empiricalQR.2d.v2(r.fitting=r.fitting,w.fitting=w.fitting,w.eval=w.eval,
                                   tau=tau,wqg=smoothness.lvl)$r.tau.wpts
    }

    # k-th score accross bins
    score.K = checkfn(quant=tau,r.eval-r0w.eval)
    return(mean(score.K,na.rm=T))

  }),na.rm=T)

  return(score)

}

quant.score.2d.checkfn.hR = function(smoothness.lvl,w,r,tau=0.95,k.folds=5,method=c("empirical","KDE"),...){
  # Perform K-fold cross-validation on the bwr parameter
  # r,w -> vectors, defined by L1 norm
  # tau -> the quantile
  # r0w -> the tau thrshold value at w

  if(min(w)<0 | max(w)>1){
    stop("angles should be on the L1 Simplex.")
  }

  grps = sample(1:k.folds,length(r),replace=T)

  score = mean(sapply(1:k.folds,function(K){
    # split data
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K]
    r.eval = r[grps==K]
    w.eval = w[grps==K]

    # fit the quantile estimates
    if(method=="KDE"){
      r0w.eval = KDE.quant.eval.2d(wpts=w.eval,r=r.fitting,w=w.fitting,tau=tau,bww=0.05,bwr=smoothness.lvl,...)
    } else if (method=="empirical"){
      r0w.eval = empiricalQR.2d.v2(r.fitting=r.fitting,w.fitting=w.fitting,w.eval=w.eval,
                                   tau=tau,wqg=smoothness.lvl)$r.tau.wpts
    }

    # k-th score accross bins
    score.K = checkfn(quant=tau,r.eval-r0w.eval)
    return(mean(score.K,na.rm=T))

  }),na.rm=T)

  return(score)

}

quant.score.checkfn = function(smoothness.lvl,w,r,tau,k.folds=5,method="KDE"){
  # r,w -> vector and matrix, defined by L1 norm
  # tau -> the quantile
  # r0w -> the tau thrshold value at w

  grps = sample(1:k.folds,length(r),replace=T)

  score = mean(sapply(1:k.folds,function(K){
    # split data
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K,]
    r.eval = r[grps==K]
    w.eval = w[grps==K,]

    # fit the quantile estimates
    if(method=="KDE"){
      r0w.eval = KDE.quant.eval(wpts=w.eval,r=r.fitting,w=w.fitting,tau=tau,bww=smoothness.lvl,bwr=0.05)
    } else{
      stop("only supports KDE")
    }
    score.K = checkfn(quant=tau,r.eval-r0w.eval)
    return(mean(score.K,na.rm=T))
  }),na.rm=T)

  return(score)
}

#########################

# A function for evaluating r_{\tau}(w) at w

KDE.quant.eval.2d = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,ker="Gaussian"){
  # r, w              -> vectors
  # bww, bwr          -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)
  # n.mesh            -> mesh for wpts
  # ker.pdf, ker.cdf  -> kernel pdf and cdf functions

  if(is.null(bww)){
    message("searching for optimal angular bandwidth for KDE threshold...")
    bww = get_bww.2d(r=r,w=w,tau=tau,bwr=bwr)
    message(paste("Fitting threshold at angular bandwidth",bww))
  }

  r.tau.wpts = sapply(wpts,function(wpts.i){
    if(ker=="Gaussian"){
      weightsw<-Gaussian.ker.pdf(w,mean=wpts.i,sd=bww)
      pos.weights = weightsw>0
      weightsw = weightsw[pos.weights]
      ccdf<-function(rc){
        ker.vals = Gaussian.ker.cdf(rc,mean=r,sd=bwr)[pos.weights]
        num = weightsw*ker.vals
        denom = weightsw
        sum(num,na.rm=T)/sum(denom,na.rm=T)
      }
    } else if(ker=="Epanechnikov"){
      weightsw<-Epanechnikov.ker.pdf(w,mean=wpts.i,sd=bww)
      pos.weights = weightsw>0
      weightsw = weightsw[pos.weights]
      ccdf<-function(rc){
        ker.vals = Epanechnikov.ker.cdf(rc,mean=r,sd=bwr)[pos.weights]
        num = weightsw*ker.vals
        denom = weightsw
        sum(num,na.rm=T)/sum(denom,na.rm=T)
      }
    }

    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this

    is_error <- FALSE
    tryCatch({
      ur<-uniroot(dummy,interval = c(0,30))$root
    },error=function(e){
      is_error <<- TRUE
    })
    if(is_error) {
      ur=NA
    }
    return(ur)
  })

  return(r.tau.wpts)
}

KDE.quant.eval = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,
                          ker.pdf=mv.Gaussian.ker.pdf,ker.cdf=Gaussian.ker.cdf){
  require(mvtnorm)

  if(is.null(bww)){
    message("searching for optimal angular bandwidth for KDE threshold...")
    bww = get_bww(r=r,w=w,tau=tau,bwr=bwr,ker.pdf=ker.pdf,ker.cdf=ker.cdf)
    message(paste("Fitting threshold at angular bandwidth",bww))
  }

  n.dims = dim(w)[2]

  r.tau.wpts = apply(wpts,1,function(wpts.i){
    if(sum(wpts.i)>1){
      return(NA)
    } else{
      weightsw = ker.pdf(x=w,mean=wpts.i,sd=bww)
      ccdf<-function(rc){
        mean(weightsw*ker.cdf(rc,mean=r,sd=bwr))/mean(weightsw)
      }
      dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
      is_error <- FALSE
      tryCatch({
        ur<-uniroot(dummy,interval = c(0,30))$root
      },error=function(e){
        is_error <<- TRUE
      })
      if(is_error) {
        ur=NA
      }
      return(ur)
    }
  })
  return(r.tau.wpts)
}

############################################

# empirical QR, fitted on 1 dataset, evaluated on another
# Also improved the notion of bin overlap

empiricalQR.2d.v2 = function (r.fitting, w.fitting, w.eval=NULL, tau = 0.95, J = 60, wqg = 1/60,
                              quantilespacing = FALSE)
{

  if(wqg<(1/J)){
    stop("Need wqg > 1/J")
  }
  wsl <- seq(0, 1 - wqg, len = J)
  wsu <- seq(wqg, 1, len = J)

  qu.r <- NULL
  for (j in 1:J) {
    re <- r.fitting[wsl[j] < w.fitting & w.fitting < wsu[j]]
    qu.r[j] <- quantile(re, tau)
  }

  if(is.null(w.eval)){
    n.mesh=200
    w.eval=seq(0,1,length.out=n.mesh)
  }
  n <- length(w.eval)
  r0w <- numeric(n)
  for (i in 1:n) {
    indvec <- rep(NA, J)
    for (j in 1:J) {
      indvec[j] <- wsl[j] < w.eval[i] & w.eval[i] < wsu[j]
    }
    r0w[i] <- mean(qu.r[indvec], na.rm = T)
  }

  return(list(wpts = w.eval, r.tau.wpts = r0w))
}

empiricalQR.3d.v2 = function (r.fitting, w.fitting, w.eval=NULL, tau = 0.95, J = 40, wg = 1/40,
                              mesh.eval=FALSE) {
  w1.fitting = w.fitting[, 1]
  w2.fitting = w.fitting[, 2]

  if(wg<(1/J)){
    stop("Need wg > 1/J")
  }
  wsl <- seq(0, 1 - wg, len = J)
  wsu <- seq(wg, 1, len = J)

  # empt=0
  qu.r <- matrix(0, J, J)
  for (j in 1:J) {
    for (k in 1:J) {
      ind <- wsl[j] < w1.fitting & w1.fitting < wsu[j] & wsl[k] < w2.fitting &
        w2.fitting < wsu[k]
      r1 <- r.fitting[ind]
      if (sum(ind) > 10) {
        qu.r[j, k] <- quantile(r1, tau)
      }
      else {
        # empt=empt+1
        # qu.r[j, k] <- NA
        qu.r[j, k] <- quantile(r1, tau)
      }
    }
  }

  if(is.null(w.eval)){
    n.mesh=100
    w.eval=expand.grid(seq(0,1,len=n.mesh),seq(0,1,len=n.mesh))
  }
  w1.eval = w.eval[, 1]
  w2.eval = w.eval[, 2]
  n <- length(w1.eval)
  r0w <- numeric(n)
  for (i in 1:n) {
    indMat <- matrix(NA, J, J)
    for (j in 1:J) {
      for (k in 1:J) {
        indMat[j, k] <- wsl[j] < w1.eval[i] & w1.eval[i] < wsu[j] &
          wsl[k] < w2.eval[i] & w2.eval[i] < wsu[k]
      }
    }
    r0w[i] <- mean(qu.r[indMat], na.rm = T)
  }

  return(list(wpts = w.eval, r.tau.wpts = r0w))
}

empiricalQR.4d.v2 = function (r.fitting, w.fitting, w.eval=NULL, tau = 0.95, J = 20, wg = 1/20,
                              mesh.eval=FALSE)
{
  w1.fitting = w.fitting[, 1]
  w2.fitting = w.fitting[, 2]
  w3.fitting = w.fitting[, 3]

  if(wg<(1/J)){
    stop("Need wg > 1/J")
  }
  wsl <- seq(0, 1 - wg, len = J)
  wsu <- seq(wg, 1, len = J)

  qu.r <- array(0, dim = c(J, J, J))
  for (j in 1:J) {
    for (k in 1:J) {
      for (l in 1:J) {
        ind <- wsl[j] < w1.fitting & w1.fitting < wsu[j] & wsl[k] <
          w2.fitting & w2.fitting < wsu[k] & wsl[l] < w3.fitting & w3.fitting < wsu[l]
        r1 <- r.fitting[ind]
        if (sum(ind) > 10) {
          qu.r[j, k, l] <- quantile(r1, tau)
        }
        else {
          qu.r[j, k, l] <- NA
        }
      }
    }
  }


  if(is.null(w.eval)){
    n.mesh=30
    w.eval=expand.grid(seq(0,1,len=n.mesh),seq(0,1,len=n.mesh),seq(0,1,len=n.mesh))
  }
  w1.eval = w.eval[, 1]
  w2.eval = w.eval[, 2]
  w3.eval = w.eval[, 3]
  n <- length(w1.eval)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = c(J, J, J))
    for (j in 1:J) {
      for (k in 1:J) {
        for (l in 1:J) {
          indArray[j, k, l] <- wsl[j] < w1.eval[i] & w1.eval[i] <
            wsu[j] & wsl[k] < w2.eval[i] & w2.eval[i] < wsu[k] &
            wsl[l] < w3.eval[i] & w3.eval[i] < wsu[l]
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }

  return(list(wpts = w.eval, r.tau.wpts = r0w))
}

empiricalQR.5d.v2 = function (r.fitting, w.fitting, w.eval=NULL, tau = 0.95, J = 10, wg = 0.15,
                              mesh.eval=FALSE)
{
  w1.fitting = w.fitting[, 1]
  w2.fitting = w.fitting[, 2]
  w3.fitting = w.fitting[, 3]
  w4.fitting = w.fitting[, 4]

  if(wg<(1/J)){
    stop("Need wg > 1/J")
  }

  # J=3;wg=0.4

  wsl <- seq(0, 1 - wg, len = J)
  wsu <- seq(wg, 1, len = J)

  qu.r <- array(0, dim = c(J, J, J, J))
  for (j in 1:J) {
    for (k in 1:J) {
      for (l in 1:J) {
        for(m in 1:J) {
          ind <- (wsl[j] < w1.fitting) & (w1.fitting < wsu[j]) &
            (wsl[k] < w2.fitting) &( w2.fitting < wsu[k]) &
            (wsl[l] < w3.fitting) &( w3.fitting < wsu[l]) &
            (wsl[m] < w4.fitting) & (w4.fitting < wsu[m])
          # print(c(,sum(ind)))
          # print(ind)
          r1 <- r.fitting[ind]
          if (sum(ind) > 2) {
            qu.r[j, k, l, m] <- quantile(r1, tau)
          }
          else {
            qu.r[j, k, l, m] <- NA
          }
        }
      }
    }
  }


  if(is.null(w.eval)){
    n.mesh=30
    w.eval=expand.grid(seq(0,1,len=n.mesh),seq(0,1,len=n.mesh),seq(0,1,len=n.mesh))
  }
  w1.eval = w.eval[, 1]
  w2.eval = w.eval[, 2]
  w3.eval = w.eval[, 3]
  w4.eval = w.eval[, 4]

  n <- length(w1.eval)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = c(J, J, J, J))
    for (j in 1:J) {
      for (k in 1:J) {
        for (l in 1:J) {
          for (m in 1:J) {
            indArray[j, k, l, m] <-
              wsl[j] < w1.eval[i] & w1.eval[i] < wsu[j] &
              wsl[k] < w2.eval[i] & w2.eval[i] < wsu[k] &
              wsl[l] < w3.eval[i] & w3.eval[i] < wsu[l] &
              wsl[m] < w4.eval[i] & w4.eval[i] < wsu[m]
          }
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }

  return(list(r0w=r0w))
  # return(list(wpts = w.eval, r.tau.wpts = r0w))
}

