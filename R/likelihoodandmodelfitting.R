# 2-d model fitting for when w is on the simplex / unit sphere. NOT Cartesian.

nll.pwlin.2d = function(psi, r, r0w, w.adj.angles, locs, norm, marg, pen.norm , pen.adj, fW.fit, joint.fit, pos.par=TRUE, pen.const=0, fixed.shape=TRUE){

  if(!fixed.shape){
    shape=psi[1]
    psi=psi[-1]
  } else {
    shape=2
  }

  # ensure we have no negative parameters
  if(pos.par) {
    if (any(psi <= 0)){
      return(1e+11)
    }
  }

  # rate <- sapply(w, gfun, par = psi)
  if(fW.fit & joint.fit){
    # exceedance radii and angles
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    G.vol = G.vol.2d(psi,locs,marg=marg)
    llw = -log(2*G.vol) -2*log(rate)
    NLL = -sum(ll1) + sum(ll2) - sum(llw)
  } else if(!fW.fit & !joint.fit){
    # exceedance radii only
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    NLL = -sum(ll1) + sum(ll2)
  } else if(fW.fit & !joint.fit){
    # exceedance angles only
    psi = c(1,psi)  # account for redundancy
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    G.vol = G.vol.2d(psi,locs,marg=marg)
    llw = -log(2*G.vol) -2*log(rate)
    NLL = -sum(llw)
  }
  if(pen.const>0){
    # pen = abs(diff(psi))
    # pen[c(1,length(pen))] = 2*pen[c(1,length(pen))]  # these appear twice in the penalty
    locs.cart = cbind(locs,1-locs)
    scaled.locs.cart = psi*locs.cart
    coplanar.vecs = apply(scaled.locs.cart,2,diff)
    norm.vecs = cbind(coplanar.vecs[,-1],-coplanar.vecs[,-2])
    grad.g = norm.vecs / apply(scaled.locs.cart[-nrow(scaled.locs.cart),] * norm.vecs, 1,sum)  # gradient of g in each region (is constant/linear in each face)
    pen = apply(grad.g,2,diff)
    if(pen.norm=="2"){
      pen = apply(pen,1,function(vec) mean(vec^2))  # used to be sum
    } else if(pen.norm=="1"){
      pen = apply(pen,1,function(vec) mean(abs(vec)))  # used to be sum
    }

    if(pen.adj){
      pen.multiplier = sapply(1:length(pen),function(r.idx){
        w.low=locs[r.idx]
        w.up=locs[r.idx+2]
        wexc.vec = sapply(w.adj.angles,function(lst) lst$w)
        num = sum(wexc.vec>=w.low & wexc.vec<=w.up)
        if(num==0){
          return(1)
        } else {
          return(num)
        }
      })
      pen = pen/pen.multiplier
    }

    pen = mean(pen)  # used to be sum
    return(NLL+(pen.const*pen))
  } else {
    return(NLL)
  }
}

opt.pwl.2d = function(NLL, r, r0w, w, locs,
                      init.val=NULL, fixed.shape=TRUE, fW.fit=FALSE,
                      joint.fit=FALSE,method="BFGS",
                      pen.const=0,pen.adj=FALSE,marg="pos",pen.norm="2",...){
  # w.adj.angles -> output of "which.adj.angles.2d"

  # if(is.null(init.val)){
  #   init.val = rep(1,length(locs))
  # }
  w.adj.angles = which.adj.angles.2d(w,locs,norm=NULL,marg=marg)
  if(fixed.shape & !fW.fit & !joint.fit){
    # only fit the conditional radial model
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
                 pen.norm=pen.norm, pen.adj=pen.adj,
                 fW.fit=F,joint.fit=F,pen.const=pen.const,fixed.shape=fixed.shape,
                 control=list(maxit=1e6),method=method,...)
    #method="L-BFGS-B",gr=NULL,lower=rep(0.01,length(locs)),upper=rep(10,length(locs))
    opt$fW.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape=2
  } else if(fixed.shape & fW.fit & joint.fit){
    # only fit the radial model and angular model togethr
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
                 pen.norm=pen.norm, pen.adj=pen.adj,
                 fW.fit=T,joint.fit=T,pen.const=pen.const,fixed.shape=fixed.shape,
                 control=list(maxit=1e6),method=method,...)
    opt$fW.par = opt$par
    opt$init.val = init.val
    opt$shape=2
  } else if(fixed.shape & fW.fit & !joint.fit){
    # only fit the angular model together, taking into account the redundancy
    init.val = init.val[-1]  # remove redundancy, fix the first parameter to 1.
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
                 pen.norm=pen.norm, pen.adj=pen.adj,
                 fW.fit=T,joint.fit=F,pen.const=pen.const,fixed.shape=fixed.shape,
                 control=list(maxit=1e6),method=method,...)
    opt$fW.par = opt$par
    opt$fW.par = c(1,opt$fW.par)  # account for redundancy
    opt$par = NULL
    opt$init.val = c(1,init.val)
    opt$shape=2
  } else if(!fixed.shape){
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
                 pen.norm=pen.norm, pen.adj=pen.adj,
                 fW.fit=F,joint.fit=F,pen.const=pen.const,fixed.shape=fixed.shape,
                 control=list(maxit=1e6),method=method,...)
    #method="L-BFGS-B",gr=NULL,lower=rep(0.01,length(locs)),upper=rep(10,length(locs))
    opt$fW.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
  }
  opt$aic = 2*(opt$value + length(opt$par))
  # return(list(mle = opt$par,fW.mle = opt$fW.par, shape=opt$shape, nll = opt$value, convergence = opt$conv,
  #             aic = opt$aic, init.val = opt$init.val))
  return(list(mle = opt$par, fW.mle = opt$fW.par, shape=2,
              nllh = opt$value, convergence = opt$conv,
              aic = opt$aic, init.val = opt$init.val))
}

# fit.pwlin.2d = function(r,w,r0w,locs,norm=NULL,marg="pos",pen.norm="2", init.val=NULL, pen.adj=F, fixed.shape=TRUE,fW.fit=FALSE,joint.fit=FALSE,pen.const=0,method="BFGS",...){
fit.pwlin.2d = function(r,w,r0w,locs,
                     init.val=NULL,fixed.shape=TRUE,fW.fit=FALSE,joint.fit=FALSE,
                     method="BFGS",pen.const=NULL,bound.fit=FALSE,
                     marg="pos",pen.norm="2",...){

  if(is.null(init.val)){
    init.val = rep(1,length(locs))
  }

  if(is.null(pen.const)){
    message("searching for gradient penalty constant...")
    pen.const = get_pen_const_2d(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,fW.fit=fW.fit,joint.fit=joint.fit)
    message(paste("Fitting with gradient penalty strength",pen.const))
  }

  gfun=function(x,par) {gfun.2d(x,par,ref.angles=locs)}

  t1 = Sys.time()

  nnodes=length(locs)
  w.adj.angles = which.adj.angles.2d(angles=w, locs, norm, marg)
  opt = opt.pwl.2d(NLL=nll.pwlin.2d, r=r, r0w=r0w,w=w,locs=locs,init.val=init.val,fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const)

  if(bound.fit){
    if(fW.fit & !joint.fit){
      stop("No need to bound anglar parameters if not jointly fitting.")
    }

    mle=opt$mle

    lik.custom = function(psi, r, r0w, w.adj.angles,
                          locs, norm, marg,
                          pen.norm, pen.adj,
                          fixed.pars, fixed.pars.idx,
                          fW.fit,joint.fit,pen.const,fixed.shape=T){

      par.full = numeric(length(locs))
      par.full[fixed.pars.idx] = fixed.pars
      par.full[-fixed.pars.idx] = psi

      return(nll.pwlin.2d(psi=par.full,r=r,r0w=r0w,w.adj.angles=w.adj.angles,
                          locs=locs,norm=norm,marg=marg,
                          pen.norm=pen.norm, pen.adj=pen.adj,
                          fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const))
    }



    # stop("generalise this to use the 1 or 2 norm")

    fixed.pars.idx = c()

    if(marg=="pos"){
      keep.fitting = any(round(as.numeric(apply(cbind(locs,1-locs)/sapply(locs,gfun,par=mle),2,max)),6)!=1)
    } else if(marg=="Rd"){
      locs.cart = pol2cart.L1Rd(w=locs)
      keep.fitting = any(round(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,min)),6)!=-1) &
        any(round(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,max)),6)!=1)
    }

    max.iters=100
    iters=1
    while(keep.fitting){
      if(marg=="pos"){
        if(0.5 %in% locs){

          idx.x  = locs >= 0.5
          idx.y  = locs <= 0.5
          idx.xy = locs == 0.5
          cm.x = max(locs/sapply(locs,gfun,par=mle))
          cm.y = max((1-locs)/sapply(locs,gfun,par=mle))
          cm.xy = 0.5/sapply(0.5,gfun,par=mle)
          if(cm.x==cm.xy & cm.y==cm.xy){
            which.cm = which(locs==0.5)
            mle[which.cm] = mle[which.cm]/cm.xy  # re-scale
            if(is.null(fixed.pars.idx)){
              fixed.pars.idx = which.cm
            } else {
              fixed.pars.idx = sort(c(fixed.pars.idx,which.cm))
            }

          } else {
            # max isn't on the diagonal
            # find the componentiwse max, and re-scale the corresponding MLEs
            idx.x = locs > 0.5
            idx.y = locs < 0.5
            unitg.x = locs/sapply(locs,gfun,par=mle)
            unitg.x[!idx.x] = 0
            unitg.y = (1-locs)/sapply(locs,gfun,par=mle)
            unitg.y[!idx.y] = 0
            which.cm.x = which.max(unitg.x)
            which.cm.y = which.max(unitg.y)
            mle[which.cm.x] = mle[which.cm.x]/(unitg.x[which.cm.x])
            mle[which.cm.y] = mle[which.cm.y]/(unitg.y[which.cm.y])

            # which mle's are re-scaled
            if(is.null(fixed.pars.idx)){
              fixed.pars.idx = sort(c(which.cm.x,which.cm.y))
            } else {
              fixed.pars.idx = sort(c(fixed.pars.idx,which.cm.x,which.cm.y))
            }
          }
        } else {
          # find the componentiwse max, and re-scale the corresponding MLEs
          idx.x = locs > 0.5
          idx.y = locs < 0.5
          unitg.x = locs/sapply(locs,gfun,par=mle)
          unitg.x[!idx.x] = 0
          unitg.y = (1-locs)/sapply(locs,gfun,par=mle)
          unitg.y[!idx.y] = 0
          which.cm.x = which.max(unitg.x)
          which.cm.y = which.max(unitg.y)
          mle[which.cm.x] = mle[which.cm.x]/(unitg.x[which.cm.x])
          mle[which.cm.y] = mle[which.cm.y]/(unitg.y[which.cm.y])

          # which mle's are re-scaled
          if(is.null(fixed.pars.idx)){
            fixed.pars.idx = sort(c(which.cm.x,which.cm.y))
          } else {
            fixed.pars.idx = unique(sort(c(fixed.pars.idx,which.cm.x,which.cm.y)))
          }
        }
      } else if(marg=="Rd"){
        unitg.coords = locs.cart/sapply(locs,gfun,par=mle)
        max.x.idx=which.max(unitg.coords[,1])
        min.x.idx=which.min(unitg.coords[,1])
        max.y.idx=which.max(unitg.coords[,2])
        min.y.idx=which.min(unitg.coords[,2])

        mle[max.x.idx] = abs(mle[max.x.idx]/(unitg.coords[max.x.idx,1]))
        mle[max.y.idx] = abs(mle[max.y.idx]/(unitg.coords[max.y.idx,2]))

        mle[min.x.idx] = abs(mle[min.x.idx]/(unitg.coords[min.x.idx,1]))
        mle[min.y.idx] = abs(mle[min.y.idx]/(unitg.coords[min.y.idx,2]))

        if(is.null(fixed.pars.idx)){
          fixed.pars.idx = sort(c(max.x.idx,max.y.idx,min.x.idx,min.y.idx))
        } else {
          fixed.pars.idx = unique(sort(c(fixed.pars.idx,max.x.idx,max.y.idx,min.x.idx,min.y.idx)))
        }
      }

      # re-fit the model
      if(fixed.shape & !fW.fit & !joint.fit){
        # only fit the conditional radial model
        init.vals=init.val[-fixed.pars.idx]  #rep(1,length(locs)-length(fixed.pars.idx))
        # opt <- optim(lik.custom, par = init.vals, r = r,
        #              r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
        #              pen.norm=pen.norm, pen.adj=pen.adj,
        #              fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,
        #              fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
        #              control=list(maxit=1e6),method=method)  #method="L-BFGS-B",gr=NULL,lower=rep(0.01,length(locs)-length(fixed.pars.idx)),upper=rep(10,length(locs)-length(fixed.pars.idx))
        opt2 = opt.pwl.2d(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                          fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                          fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx)
        opt2$fW.par = NULL  # didn't model angles
      } else if(fixed.shape & fW.fit & joint.fit){
        # only fit the radial model and angular model together
        init.vals=init.val[-fixed.pars.idx]  #rep(1,length(locs)-length(fixed.pars.idx))
        # opt <- optim(lik.custom, par = init.vals, r = r,
        #              r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, norm=norm, marg=marg,
        #              pen.norm=pen.norm, pen.adj=pen.adj,
        #              fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,
        #              fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
        #              control=list(maxit=1e6),method=method)
        opt2 = opt.pwl.2d(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                          fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                          fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx)
        opt2$fW.par = opt2$mle
      } else {
        stop("ERROR in likelihood setup.")
      }

      mle[-fixed.pars.idx] = opt2$mle

      if(marg=="pos"){
        locs.cart = cbind(locs,1-locs)
        keep.fitting = any(round(as.numeric(apply(cbind(locs,1-locs)/sapply(locs,gfun,par=mle),2,max)),6)!=1)
      } else if(marg=="Rd"){
        locs.cart = pol2cart.L1Rd(w=locs)
        keep.fitting = any(round(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,min)),6)!=-1) &
          any(round(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,max)),6)!=1)
      }

      if(!keep.fitting){
        # message(paste("comp. min:",paste(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,min)),collapse=" ")))
        # message(paste("comp. max:",paste(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,max)),collapse=" ")))
        break
      }

      iters=iters+1
      if(iters>max.iters){
        message("Max iters in bounding algorithm reached.")
        break
      }

    }
    # return(list(mle=mle,
    #             idx=fixed.pars.idx,
    #             nll=opt$value))
    t2 = Sys.time()

    if(joint.fit){
      fW.mle = mle
    } else {
      fW.mle = NULL
    }

    return(list(mle =mle,
                fixed.pars.idx=fixed.pars.idx,
                fW.mle = fW.mle,
                shape=2,
                pen.const=pen.const,
                nll = NULL,
                convergence = opt2$conv,
                init.val = init.val,
                info=t2-t1))
  } else {
    t2 = Sys.time()

    return(list(mle =opt$mle,
                fixed.pars.idx=NULL,
                fW.mle = opt$fW.mle,
                shape=2,
                pen.const=pen.const,
                nll = opt$nllh,
                convergence = opt$conv,
                init.val = init.val,
                info=t2-t1))
  }
}

get_pen_const_2d = function(r,w,r0w,locs,init.val,fW.fit=FALSE,joint.fit=FALSE){
  if(fW.fit & !joint.fit){
    pen.consts = seq(0,30,length.out=100)
  } else {
    pen.consts = seq(0,4,length.out=100)
  }

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  pen.consts.scores = rep(NA,length(pen.consts))

  iter=1
  for(pen.const.idx in c(1:length(pen.consts))){

    pen.const = pen.consts[pen.const.idx]
    # pen.const.W = pen.consts.W[pen.const.idx]

    # lik.scores.R13 = rep(NA,k.folds)

    # for(K in 1:k.folds){
    K = 1
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K]
    r0w.fitting = r0w[grps!=K]
    r.eval = r[grps==K]
    w.eval = w[grps==K]
    r0w.eval = r0w[grps==K]

    is_error <- FALSE
    tryCatch({
      mod.fit = fit.pwlin.2d(r=r.fitting,r0w=r0w.fitting,w=w.fitting,
                               locs=locs,pen.const=pen.const,method="BFGS",
                               init.val=init.val,
                             fW.fit=fW.fit,joint.fit=joint.fit,bound.fit=FALSE)
    },error=function(e){
      is_error <<- TRUE
    })
    if(is_error) {
      # print(pen.const)
      pen.consts.scores[pen.const.idx] = 1e10
      next
    }

    w.eval.adj.angles = which.adj.angles.2d(angles=w.eval,locs=locs)

    if(fW.fit & !joint.fit){
      psi.vals = mod.fit$fW.mle
    } else {
      psi.vals = mod.fit$mle
    }

    pen.consts.scores[pen.const.idx] =
      nll.pwlin.2d(psi=psi.vals, r=r.eval, r0w=r0w.eval,
                                     w.adj.angles=w.eval.adj.angles,
                                     locs=locs,
                                     norm="1", marg="pos",pen.norm=NULL, pen.adj=NULL,
                                     fW.fit=fW.fit,joint.fit=joint.fit)
  }
  return(pen.consts[which.min(pen.consts.scores)])
}

######################################################################

# general d-dimensional model fitting for when w is on the simplex.

# pwl.L1.pen = function(psi,locs){
#   # Experimental L1 penalty
#
#   require(geometry)
#   pen.vals = sapply(1:length(psi), function(i){
#     # print(locs[i,])
#     psi.i = psi[i]
#     locs.i = locs[i,]
#     psi.rest = psi[-i]
#     locs.rest = locs[-i,]
#     suppressWarnings({
#       g.val.i = gfun.pwl(x=locs.i,par = psi.rest,ref.angles=locs.rest)
#     })
#     return(abs(psi.i - (1/g.val.i)))
#   })
#   return(mean(pen.vals,na.rm=T))  # used to be sum
# }
pwl.L2.pen = function(psi,locs,del.tri){
  require(geometry)
  tri = del.tri$tri

  locs.scaled = locs*psi

  num.cols = dim(locs)[2]

  grad.g = lapply(1:nrow(tri),function(row.idx){
    idx = tri[row.idx,]
    idx.fix = idx[1]
    idx.rest = idx[-1]
    coplanar.mat = do.call(rbind,lapply(idx.rest, function(ii){
      return(locs.scaled[idx.fix,] - locs.scaled[ii,])
    }))
    # norm.vec = c(det(coplanar.mat[,-1]),-det(coplanar.mat[,-2]),det(coplanar.mat[,-3]),-det(coplanar.mat[,-4]))
    norm.vec = suppressWarnings({
      c(1,-1)*sapply(c(1:num.cols),function(i){det(coplanar.mat[,-i])})
    })
    grad.g.val = norm.vec / sum(norm.vec * locs.scaled[idx.fix,])
    return(grad.g.val)
  })
  # pen = 0
  idx.vecs = NULL
  for(i in 1:nrow(locs)){
    which.has.i = apply(tri,1,function(ii){
      i %in% ii
    })
    if(sum(which.has.i)==1){
      next  # only 1 face, can't compare gradient to another
    } else if(sum(which.has.i)>1){
      which.grad.gs = (1:length(grad.g))[which.has.i]
      idx.combs = combn(1:length(which.grad.gs),2)
      idx.vec = do.call(rbind,lapply(1:ncol(idx.combs),function(i.col){
        return(sort(which.grad.gs[idx.combs[,i.col]]))
      }))
      idx.vecs = rbind(idx.vecs,idx.vec)
    }
  }
  idx.vecs = idx.vecs[!duplicated(idx.vecs), ]
  grad.g.mat = do.call(rbind,grad.g)
  pen = apply(idx.vecs,1,function(ii){
    pp = diff(grad.g.mat[ii,])
    return(mean(pp^2))  # used to be sum
  })
  return(mean(pen))   # used to be sum
}

nll.pwlin = function(psi, r, r0w, w.adj.angles, del.tri, locs, fW.fit, joint.fit, pen.const=0, pos.par=TRUE){

  # ensure we have no negative parameters
  if(pos.par) {
    if (any(psi <= 0)){
      return(1e+11)
    }
  }

  # # REMOVE THIS...
  # if (any(psi > 5)){
  #   return(1e+11)
  # }

  num.cols = length(w.adj.angles[[1]]$w)

  if(fW.fit & joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = num.cols, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = num.cols, rate = rate, lower.tail = F))
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(ll1) + sum(ll2) - sum(llw)
  } else if(!fW.fit & !joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = num.cols, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = num.cols, rate = rate, lower.tail = F))
    NLL = -sum(ll1) + sum(ll2)
  } else if(fW.fit & !joint.fit){
    psi = c(1,psi)  # account for redundancy
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(llw)
  }
  # if(pen.const.L1>0){
  #   L1.pen.val = pwl.L1.pen(psi=psi,locs=locs)
  #   NLL = NLL+(pen.const.L1*L1.pen.val)
  # }
  if(pen.const>0){
    L2.pen.val = pwl.L2.pen(psi=psi,locs=locs,del.tri=del.tri)
    NLL = NLL+(pen.const*L2.pen.val)
  }
  # print(psi)
  # print(NLL)
  return(NLL)
}

opt.pwl = function(NLL, r, r0w, w, locs,
                     init.val=NULL, fixed.shape=TRUE, fW.fit=FALSE,
                     joint.fit=FALSE,method="BFGS",
                     pen.const=0,...){

  # t1 = Sys.time()

  if(is.null(init.val)){
    init.val = rep(1,nrow(locs))
  }

  w.adj.angles = which.adj.angles(w,locs)

  num.cols = length(w.adj.angles[[1]]$w)
  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)

  if(fixed.shape & !fW.fit & !joint.fit){
    # fit the conditional radial model only
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 fW.fit=F,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = NULL  # didn't model angles
  } else if(fixed.shape & fW.fit & joint.fit){
    # fit the radial model and angular model together
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 fW.fit=T,joint.fit=T,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = opt$par
  } else if(fixed.shape & fW.fit & !joint.fit){
    # fit the angular model only
    init.val = init.val[-1]  # remove redundancy, fit the first parameter to 1.
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 fW.fit=T,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = opt$par
    opt$fW.par = c(1,opt$fW.par)  # account for redundancy
    opt$par = NULL
  } else if(!fixed.shape){
    stop("fitting shape not implemented.")
  }

  # t2 = Sys.time()

  opt$init.val = init.val#c(3,init.val)
  opt$aic = 2*(opt$value + length(opt$par))
  return(list(mle = opt$par, fW.mle = opt$fW.par, shape=dim(w)[2],
              nllh = opt$value, convergence = opt$conv,
              aic = opt$aic, init.val = opt$init.val))
}

fit.pwlin = function(r,w,r0w,locs,
                     init.val=NULL,fixed.shape=TRUE,fW.fit=FALSE,joint.fit=FALSE,
                     method="BFGS",pen.const=NULL,bound.fit=FALSE,...){

  if(fW.fit & !joint.fit & bound.fit){
    stop("No need to bound anglar parameters if not jointly fitting.")
  }

  t1 = Sys.time()

  if(is.null(init.val)){
    init.val = rep(1,nrow(locs))
  }

  if(is.null(pen.const)){
    message("searching for gradient penalty constant...")
    pen.const = get_pen_const(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,fW.fit=fW.fit,joint.fit=joint.fit)
    message(paste("Fitting with gradient penalty strength",pen.const))
  }

  num.cols = dim(w)[2]

  w.adj.angles = which.adj.angles(angles=w,locs=locs)
  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)

  # fixed.pars.idx=NULL
  opt <- opt.pwl(NLL=nll.pwlin,r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,...)

  # print(opt)

  if(bound.fit){

    mle = opt$mle

    lik.custom = function(psi, r, r0w, w.adj.angles, locs = locs, fixed.pars,
                          fixed.pars.idx, fW.fit, joint.fit, pos.par = TRUE,
                          pen.const=0,del.tri=del.tri){

      psi.full = numeric(nrow(locs))
      psi.full[fixed.pars.idx] = fixed.pars
      psi.full[-fixed.pars.idx] = psi

      # if(pos.par) {
      #   if (any(psi.full <= 0)){
      #     return(1e+11)
      #   }
      # }

      return(nll.pwlin(psi=psi.full,r=r,r0w=r0w,w.adj.angles=w.adj.angles,
                       locs=locs, fW.fit=fW.fit, joint.fit=joint.fit, pos.par=TRUE,
                       pen.const=pen.const,del.tri=del.tri))

    }
    fixed.pars.idx = c()
    is.bounded=F   # kind of like a dummy variable, will terminate the loop inside
    while(!is.bounded){
      # print(i)

      # evaluate at nodes
      unitg = locs/gfun.pwl(x=locs,par=mle,ref.angles=locs)#apply(locs,1,gfun.4d,par=mle,ref.angles=locs)  # need to round here

      # what is max along nodes
      max.vals = apply(unitg,2, max)
      which.max.vals = apply(unitg,2, which.max)  # row values

      max.df = data.frame(coord=c(1:num.cols),max=max.vals,loc=which.max.vals)
      max.lst = lapply(1:nrow(locs),function(loc){
        cond = max.df$loc==loc
        if(all(cond==F)){
          return(NULL)
        } else {
          return(list(loc=loc,
                      val=max(max.df$max[cond])))
        }
      })
      max.lst = max.lst[lapply(max.lst,length)>0]

      if(all(round(max.vals,num.cols)==1)){
        break  # algorithm complete, all comp. max. vals are 1
      }
      for(lst in max.lst){
        fixed.pars.idx = unique(sort(c(fixed.pars.idx,lst$loc)))
        mle[lst$loc] = mle[lst$loc]/lst$val
      }

      if(nrow(locs)==length(fixed.pars.idx)){
        break  # we fixed all the parameters
      }

      if(fixed.shape & !fW.fit & !joint.fit){
        # fit the conditional radial model only
        init.vals = init.val[-fixed.pars.idx]
        # opt <- optim(lik.custom, par = init.vals, r = r,
        #              r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
        #              fW.fit=F,joint.fit=F,pen.const=pen.const,del.tri=del.tri,
        #              fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
        #              control=list(maxit=1e6,reltol=1e-5),method=method,...)
        opt2 = opt.pwl(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                      fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                      fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,...)
        opt2$fW.par = NULL  # didn't model angles
      } else if(fixed.shape & fW.fit & joint.fit){
        # fit the radial model and angular model together
        init.vals = init.val[-fixed.pars.idx]
        # opt <- optim(lik.custom, par = init.vals, r = r,
        #              r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
        #              fW.fit=T,joint.fit=T,pen.const=pen.const,del.tri=del.tri,
        #              fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
        #              control=list(maxit=1e6,reltol=1e-5),method=method,...)
        opt2 = opt.pwl(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                       fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                       fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,...)
        opt2$fW.par = opt2$mle
      } else {
        stop("ERROR in likelihood setup")
      }
      mle[-fixed.pars.idx] = opt2$mle
      # i=i+1
    }

    t2 = Sys.time()

    if(joint.fit){
      fW.mle = mle
    } else {
      fW.mle = NULL
    }

    return(list(mle = mle, fW.mle = fW.mle, shape=dim(w)[2],
                pen.const=pen.const,
                fixed.pars.idx=fixed.pars.idx,
                nllh = NULL, convergence = opt2$conv,
                aic = NULL, init.val = init.val,
                info = t2-t1))
  } else {
    t2 = Sys.time()
    opt$pen.const=pen.const
    opt$info = t2-t1
    opt$fixed.pars.idx=NULL
    return(opt)
  }

  # return(list(mle=mle,
  #             fixed.pars.idx=fixed.pars.idx))

  # t2 = Sys.time()

  # return(list(mle = mle, fW.mle = opt$fW.par, shape=dim(w)[2],
  #             nllh = opt$value, convergence = opt$conv,
  #             aic = opt$aic, init.val = opt$init.val,
  #             info = t2-t1))
}


get_pen_const = function(r,w,r0w,locs,init.val,fW.fit=FALSE,joint.fit=FALSE){
  if(fW.fit & !joint.fit){
    pen.consts = seq(0,30,length.out=30)
  } else {
    pen.consts = seq(0,4,length.out=30)
  }

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  pen.consts.scores = rep(NA,length(pen.consts))

  iter=1
  for(pen.const.idx in c(1:length(pen.consts))){
    # print(pen.consts.scores)


    pen.const = pen.consts[pen.const.idx]
    # pen.const.W = pen.consts.W[pen.const.idx]

    # lik.scores.R13 = rep(NA,k.folds)

    # for(K in 1:k.folds){
    K = 1
    r.fitting = r[grps!=K]
    w.fitting = w[grps!=K,]
    r0w.fitting = r0w[grps!=K]
    r.eval = r[grps==K]
    w.eval = w[grps==K,]
    r0w.eval = r0w[grps==K]

    is_error <- FALSE
    tryCatch({
      mod.fit = fit.pwlin(r=r.fitting,r0w=r0w.fitting,w=w.fitting,
                             locs=locs,pen.const=pen.const,method="BFGS",
                             init.val=init.val,
                             fW.fit=fW.fit,joint.fit=joint.fit,bound.fit=FALSE)
    },error=function(e){
      is_error <<- TRUE
    })
    if(is_error) {
      # print(pen.const)
      pen.consts.scores[pen.const.idx] = 1e10
      next
    }

    w.eval.adj.angles = which.adj.angles(angles=w.eval,locs=locs)

    if(fW.fit & !joint.fit){
      psi.vals = mod.fit$fW.mle
    } else {
      psi.vals = mod.fit$mle
    }

    pen.consts.scores[pen.const.idx] =
      nll.pwlin(psi=psi.vals, r=r.eval, r0w=r0w.eval,
                   w.adj.angles=w.eval.adj.angles,
                   locs=locs,
                   fW.fit=fW.fit,joint.fit=joint.fit)
  }
  return(pen.consts[which.min(pen.consts.scores)])
}

