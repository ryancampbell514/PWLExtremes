pwl.L2.pen.trianglewise = function(psi,locs,del.tri,adj.list){
  require(geometry)
  tri = del.tri$tri

  num.tris = nrow(tri)

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

  grad.g.mat = do.call(rbind,grad.g)

  pen = sapply(1:length(adj.list), function(k){
    adj.list.k = adj.list[[k]]
    pen.k = lapply(adj.list.k, function(j){
      diff(grad.g.mat[c(j,k),])
    })
    pen.k = mean(apply(do.call(rbind,pen.k)^2,1,sum))
    return(pen.k)
  })
  return(mean(pen))   # should this be "mean" or "sum"?
}

nll.pwlin.TEST = function(psi, r, r0w, w.adj.angles, fixed.shape, del.tri, locs,
                          adj.list,fW.fit, joint.fit, pen.const=0, pos.par=TRUE){

  num.cols = length(w.adj.angles[[1]]$w)

  if(!fixed.shape){
    shape=psi[1]
    psi=psi[-1]
    if(shape < 0){
      return(1e+11)
    }
  } else {
    shape=num.cols
  }

  # ensure we have no negative parameters
  if(pos.par) {
    if (any(psi <= 0)){
      return(1e+11)
    }
  }

  if(fW.fit & joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(ll1) + sum(ll2) - sum(llw)
  } else if(!fW.fit & !joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    NLL = -sum(ll1) + sum(ll2)
  } else if(fW.fit & !joint.fit){
    psi = c(1,psi)  # account for redundancy
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(llw)
  }
  if(pen.const>0){
    L2.pen.val = pwl.L2.pen.trianglewise(psi=psi,locs=locs,del.tri=del.tri,adj.list=adj.list)
    NLL = NLL+(pen.const*L2.pen.val)
  }
  print(psi)
  print(NLL)
  return(NLL)
}

opt.pwl.TEST = function(NLL, r, r0w, w, locs, w.adj.angles, del.tri, adj.list,
                        init.val=NULL, fixed.shape=TRUE, fW.fit=FALSE,
                        joint.fit=FALSE,method="BFGS",
                        pen.const=0,...){

  # w.adj.angles = which.adj.angles(w,locs)

  num.cols = length(w.adj.angles[[1]]$w)
  # del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  # adj.list = ij.couples(locs)

  if(fixed.shape & !fW.fit & !joint.fit){
    # fit the conditional radial model only
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 adj.list=adj.list,
                 fW.fit=F,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixed.shape=fixed.shape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape=num.cols
  } else if(fixed.shape & fW.fit & joint.fit){
    # fit the radial model and angular model together
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 adj.list=adj.list,
                 fW.fit=T,joint.fit=T,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixed.shape=fixed.shape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = opt$par
    opt$init.val = init.val
    opt$shape=num.cols
  } else if(fW.fit & !joint.fit){
    # fit the angular model only
    init.val = init.val[-1]  # remove redundancy, fit the first parameter to 1.
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 adj.list=adj.list,
                 fW.fit=T,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixed.shape=fixed.shape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$fW.par = opt$par
    opt$fW.par = c(1,opt$fW.par)  # account for redundancy
    opt$par = NULL
    opt$init.val = c(1,init.val)
    opt$shape=NULL
  } else if(!fixed.shape & !fW.fit & !joint.fit){
    # fit the conditional radial model only
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 adj.list=adj.list,
                 fW.fit=F,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixed.shape=fixed.shape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
  } else if(!fixed.shape & fW.fit & joint.fit){
    # fit the radial model and angular model together
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 adj.list=adj.list,
                 fW.fit=T,joint.fit=T,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixed.shape=fixed.shape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,...)
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
    opt$fW.par = opt$par
  }

  opt$init.val = init.val#c(3,init.val)
  opt$aic = 2*(opt$value + length(opt$par))
  return(list(mle = opt$par, fW.mle = opt$fW.par, shape=opt$shape,
              nllh = opt$value, convergence = opt$conv,
              aic = opt$aic, init.val = opt$init.val))
}

fit.pwlin.TEST = function(r,w,r0w,locs,
                     init.val=NULL,fixed.shape=TRUE,fW.fit=FALSE,joint.fit=FALSE,
                     method="BFGS",pen.const=NULL,bound.fit=FALSE,...){

  if(fW.fit & !joint.fit & bound.fit){
    stop("No need to bound anglar parameters if not jointly fitting.")
  }

  t1 = Sys.time()

  # if(is.null(init.val)){
  #   init.val = rep(1,nrow(locs))
  # }
  if(is.null(init.val) & fixed.shape){
    init.val = rep(1,nrow(locs))
  } else if(is.null(init.val) & !fixed.shape){
    init.val = rep(1,nrow(locs)+1)
  }

  if(is.null(pen.const)){
    message("searching for gradient penalty constant...")
    pen.const = get_pen_const(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,fW.fit=fW.fit,joint.fit=joint.fit)
    message(paste("Fitting with gradient penalty strength",pen.const))
  }

  num.cols = dim(w)[2]

  w.adj.angles = which.adj.angles(angles=w,locs=locs)
  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  adj.list = adj.DT(locs)

  # fixed.pars.idx=NULL
  opt <- opt.pwl.TEST(NLL=nll.pwlin.TEST,r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,
                 w.adj.angles=w.adj.angles, del.tri=del.tri, adj.list=adj.list,
                 fixed.shape=fixed.shape,fW.fit=fW.fit,joint.fit=joint.fit,
                 pen.const=pen.const,method=method,...)

  # print(opt)

  if(bound.fit){

    shape.val = opt$shape
    mle = opt$mle

    lik.custom = function(psi, r, r0w, w.adj.angles, locs = locs, adj.list=adj.list,
                          fixed.pars,fixed.pars.idx, fixed.shape, fW.fit, joint.fit, pos.par = TRUE,
                          pen.const=0,del.tri=del.tri){

      # psi.full = numeric(nrow(locs))
      # psi.full[fixed.pars.idx] = fixed.pars
      # psi.full[-fixed.pars.idx] = psi
      if(fixed.shape){
        par.full = numeric(nrow(locs))
        par.full[fixed.pars.idx] = fixed.pars
        par.full[-fixed.pars.idx] = psi
      } else {
        par.full = numeric(nrow(locs)+1)
        par.full[fixed.pars.idx+1] = fixed.pars
        par.full[-(fixed.pars.idx+1)] = psi
      }
      # if(pos.par) {
      #   if (any(psi.full <= 0)){
      #     return(1e+11)
      #   }
      # }

      return(nll.pwlin.TEST(psi=par.full,r=r,r0w=r0w,w.adj.angles=w.adj.angles,
                       locs=locs, adj.list=adj.list, fW.fit=fW.fit,
                       joint.fit=joint.fit, pos.par=TRUE,
                       pen.const=pen.const,del.tri=del.tri,fixed.shape=fixed.shape))

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

      init.vals=mle[-fixed.pars.idx]
      if(!fixed.shape){
        init.vals = c(shape.val,init.vals)
      }

      if(nrow(locs)==length(fixed.pars.idx)){
        break  # we fixed all the parameters
      }

      if(!fW.fit & !joint.fit){
        # fit the conditional radial model only
        init.vals = init.val[-fixed.pars.idx]
        opt2 = opt.pwl.TEST(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                       w.adj.angles=w.adj.angles, del.tri=del.tri, adj.list=adj.list,
                       fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                       fixed.shape=fixed.shape,
                       fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,...)
        opt2$fW.par = NULL  # didn't model angles
      } else if(fW.fit & joint.fit){
        # fit the radial model and angular model together
        init.vals = init.val[-fixed.pars.idx]
        opt2 = opt.pwl.TEST(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                       w.adj.angles=w.adj.angles, del.tri=del.tri, adj.list=adj.list,
                       fW.fit=fW.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                       fixed.shape=fixed.shape,
                       fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,...)
        opt2$fW.par = opt2$mle
      } else {
        stop("ERROR in likelihood setup")
      }
      mle[-fixed.pars.idx] = opt2$mle
      shape.val = opt2$shape
      # i=i+1
    }

    t2 = Sys.time()

    if(joint.fit){
      fW.mle = mle
    } else {
      fW.mle = NULL
    }

    return(list(mle = mle, fW.mle = fW.mle, shape=shape.val,
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
}
