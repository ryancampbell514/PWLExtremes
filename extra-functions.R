# Functions that are not yet in the geometricMVE package but are required
# for reproducing certain results in the piecewise-linear paper.

# Threshold projection:
# call this *.RC because I need it, but it is not exported in geometricMVE
KDE.thresh.eval.RC = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){
  require(mvtnorm)
  
  n.dims = dim(w)[2]
  
  r.tau.wpts = apply(wpts,1,function(wpts.i){
    if(sum(wpts.i)>1){
      return(NA)
    } else{
      if(length(bww)==1){
        weightsw = dmvnorm(w[,-n.dims],mean=wpts.i[-n.dims],sigma=(bww^2)*diag(n.dims-1))
      } else {
        weightsw.i = lapply(1:ncol(bww),function(bww.idx){
          dnorm(x = w[,bww.idx],mean=wpts.i[bww.idx], sd = bww[,bww.idx])
        })
        weightsw = apply(do.call(cbind,weightsw.i),1,prod)
      }
      
      ccdf<-function(rc){
        mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
      }
      dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
      is_error <- FALSE
      tryCatch({
        ur<-uniroot(dummy,interval = c(0,up))$root
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

#' Plot fitted threshold for d>3

#' @param thresh.fit output from fit.thresh
#' @param resolution number giving plotting resolution of fitted threshold.
#' @param add logical - add to existing plot?
#' @param which.proj a vector of integers of length d-3 indicating which variables to project
#' @return rgl plot of fitted threshold over data sample
#'
#' @export
#' 
plotfittedthresh.3dproj=function(thresh.fit, resolution=30, add=FALSE, which.proj,
                                 xlab="x1",ylab="x2",zlab="x3"){
  d = dim(thresh.fit$w)[2]  # thresh.fit$w has the last column removed
  if(d-length(which.proj)!=3){
    stop(paste("which.proj needs to be a vector of length",d-3))
  }
  
  rtilde.proj<-function(x,which.proj,upper=30,res=100)
  {
    if(any(is.na(x))){
      return(NA)
    }
    
    dummy<-function(y){
      w.inp = rep(NA,d)
      w.inp[which.proj] = y
      w.inp[-which.proj] = x
      
      # thresh.fit2 = thresh.fit
      # thresh.fit2$w = cbind(thresh.fit2$w,1-apply(thresh.fit2$w,1,sum))
      
      # print(matrix(w.inp/sum(w.inp),nrow=1))
      # stop()
      
      sum(w.inp) / KDE.thresh.eval.RC(wpts = matrix(w.inp/sum(w.inp),nrow=1),
                                      r = thresh.fit$r, w=thresh.fit$w,
                                      bww=thresh.fit$bww, tau=thresh.fit$tau)
    }
    if(d==4){
      opt<-optimize(dummy,interval=c(0,upper))
      return(opt$objective)
    } else{
      # print(x)
      opt<-optim(par=rep(1,length(which.proj)),fn=dummy)
      return(opt$value)
    }
  }
  
  # the plotting angles
  wseq<-seq(0,1,len=resolution)
  wgrid = expand.grid(wseq,wseq)
  wmat<-cbind(wgrid,1-apply(wgrid,1,sum))
  nacond = wmat[,3] <0
  wmat[nacond,3] = NA
  
  thresh.proj<-apply(wmat,1,rtilde.proj,which.proj=which.proj)
  thresh.proj[nacond] = NA
  thresh.proj.mat = matrix(thresh.proj,sqrt(length(thresh.proj)),sqrt(length(thresh.proj)))
  
  full.ds = thresh.fit$r * cbind(thresh.fit$w,1-apply(thresh.fit$w,1,sum))
  
  usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                     -0.09606559,-0.2029479,0.974465668,0,
                     0.38427672,0.8955331,0.224392071,0,
                     0,0,0,1)
                   ,4,4,byrow=T) 
  # axis.names = paste0("x",as.character(c(1:d)[-which.proj]))
  if(!add){
    open3d()
    plot3d(full.ds[,-which.proj],col="black",alpha=1,
           xlab=xlab,ylab=ylab,zlab=zlab)
  }
  surface3d(wmat[,1]/thresh.proj.mat,
            wmat[,2]/thresh.proj.mat,
            wmat[,3]/thresh.proj.mat,
            col="red",alpha=0.4)
  axes3d()
  view3d(userMatrix = usermat)
}

###########################################################################

# CODE TO IMPLEMENT THE KDE THRESHOLD ON THE SIMPLEX, WITH THE OPTION FOR 
# ADAPTIVE BANDWIDTH SELECTION.
# EXPERIMENT WITH USING THE EMPIRICAL CDF FOR THE RADIAL CDF

#' Find a high threshold of R|W which is (approximately) the tau-quantile of this variable
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on the unit simplex)
#' @param tau level at which to calculate threshold
#' @param method character string specifying either "empirical" for using a binning method, or "KDE" for the method based on kernel density estimation as outlined in Campbell and Wadsworth (2025). Defaults to "KDE".
#' @param bin.mesh numerical value affecting number of bins for estimation when method="empirical"
#' @param overlap numerical value affecting overlap of bins for estimation when method="empirical"
#' @param bww a scalar bandwidth of the angular kernel density when method="KDE". If NULL, adaptive selection if performed. (default NULL)
#' @param up numerical value giving upper limit for root finding of r0w when method="KDE" (default 30)
#' @param alpha numerical value in (0,1]  when method="KDE" (default 0.5)

#' @return list containing elements r0w (estimated threshold for each given w), r, w, tau, bww, and method
#' @export

fit.thresh = function(r,w,tau=0.95,method="KDE",bin.mesh=NULL,overlap=NULL,bww=NULL,up=30,alpha=0.5){
  if(any(w<0)){stop("values of w must be in the unit simplex")}
  if(is.vector(w)){
    w<-cbind(w,1-w)
  } else{
    sw<-apply(w,1,sum)
    sw = round(sw,4)
    if(all(sw<1)){
      w<-cbind(w,1-sw)
    }
  }
  
  if(dim(w)[2]==2 && method=="KDE"){
    # w = w[,-ncol(w)]
    require(pdfCluster)
    if(is.null(bww)){
      require(pdfCluster)
      bww = as.numeric(pdfCluster::kepdf(x = w[,-ncol(w)], kernel = "gaussian", bwtype="adaptive",
                                         alpha=alpha)@par$hx)  # can play with the apha hyperparameter
    }
    r0w = radial.thresh.KDE.2d(r=r,w=w[,-ncol(w)],tau=tau,bww=bww,up=up)
  } else if(dim(w)[2]>2 && method=="KDE") {
    # w = w[,-ncol(w)]
    if(is.null(bww)){
      require(pdfCluster)
      bww = pdfCluster::kepdf(x = w[,-ncol(w)], kernel = "gaussian", bwtype="adaptive",
                              alpha=alpha)@par$hx
    } else if (length(bww)==1) {
      bww = matrix(bww,nrow(w),ncol(w)-1)  # make the scalar value a matrix to be compatible with the code in radial.thresh.KDE
    } else {stop("bww must be NULL or a single numeric value")}
    r0w = radial.thresh.KDE(r=r,w=w[,-ncol(w)],tau=tau,bww=bww,up=up)
  } else if(dim(w)[2]==2 && method=="empirical"){
    # r0w = emp.thresh(r=r,w=w,tau=tau)
    emp.thresh.output = emp.thresh.2d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==3 && method=="empirical"){
    emp.thresh.output = emp.thresh.3d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==4 && method=="empirical"){
    emp.thresh.output = emp.thresh.4d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==5 && method=="empirical"){
    emp.thresh.output = emp.thresh.5d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]>5 && method=="empirical"){
    stop("Empirical threshold estimation is not yet implemented for d>5. Use the KDE method instead.")
  }
  if(method=="KDE"){
    return(list(r0w=r0w,
                w=w,
                r=r,
                tau=tau,
                bww=bww,
                method=method))
  } else if(method=="empirical"){
    return(list(r0w=emp.thresh.output$r0w,
                w=w,
                r=r,
                tau=tau,
                quant.grid=emp.thresh.output$quant.grid,  # needed for eval.thresh.emp...
                method=method,
                bin.mesh=bin.mesh,
                overlap=overlap))
  }
}

#' Evaluate a high threshold of R|W at a set of new angles W
#' 
#' @param fit.thresh.out the output of the 'fit.thresh' function
#' @param w.eval a vector of length n or a matrix with n rows
#' @param up numerical value giving upper limit for root finding of r0w when using KDE method (default 30)
#' 
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @export
eval.thresh = function(fit.thresh.out,w.eval,up=30){
  
  # if(any(w.eval<0)){stop("values of w must be in the unit simplex")}
  if(is.vector(w.eval)){
    w.eval<-cbind(w.eval,1-w.eval)
  } else{
    sw<-apply(w.eval,1,sum)
    sw = round(sw,4)
    if(all(sw<1)){
      w.eval<-cbind(w.eval,1-sw)
    }
  }
  # print(w.eval)
  
  if(dim(w.eval)[2]==2 & fit.thresh.out$method=="KDE"){
    w.eval = w.eval[,1]
    r0w = KDE.thresh.eval.2d(wpts=w.eval,
                             r=fit.thresh.out$r,w=fit.thresh.out$w,up=up,
                             tau=fit.thresh.out$tau,bww=fit.thresh.out$bww)
  } else if(dim(w.eval)[2]>2 & fit.thresh.out$method=="KDE") {
    w.eval = w.eval[,-ncol(w.eval)]
    r0w = KDE.thresh.eval(wpts=w.eval,
                          r=fit.thresh.out$r,w=fit.thresh.out$w,
                          up=up,tau=fit.thresh.out$tau,bww=fit.thresh.out$bww)
  } else if(dim(w.eval)[2]==2 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.2d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==3 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.3d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==4 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.4d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==5 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.5d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]>5 & fit.thresh.out$method=="empirical"){
    stop("Empirical threshold estimation is not yet implemented for d>5. Use the KDE method instead.")
  }
  
  return(r0w)
  
}

############################################################################

# KDE threshold estimation

#' d=2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using a KDE approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable
#' @param tau level at which to calculate threshold
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
radial.thresh.KDE.2d = function(r,w,tau=0.95,bww,bwr=0.05,up=30){
  
  r0w = sapply(w, function(ww){
    weightsw<-dnorm(w,mean=ww,sd=bww)
    
    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
      # weightsr = sapply(r,function(ri) mean())  # how to do empirical K_R?
      # mean(weightsw*mean(rc<r))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,up))
    return(ur$root)
  })
  
  return(r0w)
}

#' d=2 setting: Evaluate a high threshold of R|W at a set of new angles W using a KDE method
#' 
#' @param wpts a vector of length n or a matrix with n rows and 2 columns
#' @param r Values of radial variable used to create the KDE threshold
#' @param w values of angular variable used to create the KDE threshold
#' @param tau level at which to calculate threshold
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
KDE.thresh.eval.2d = function(wpts,r,w,tau=0.95,bww=NULL,bwr=0.05,up=30){
  # r, w              -> vectors
  # bwr          -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)
  # n.mesh            -> mesh for wpts
  # ker.pdf, ker.cdf  -> kernel pdf and cdf functions
  
  # require(pdfCluster)
  # if(is.null(bww)){
  #   bww = as.numeric(pdfCluster::kepdf(x = w, kernel = "gaussian", bwtype="adaptive",
  #                                      alpha=0.5)@par$hx)
  # }
  r.tau.wpts = sapply(wpts,function(wpts.i){
    weightsw<- dnorm(w,mean=wpts.i,sd=bww)
    ccdf<-function(rc){
      ker.vals = pnorm(rc,mean=r,sd=bwr)
      num = weightsw*ker.vals
      denom = weightsw
      sum(num,na.rm=T)/sum(denom,na.rm=T)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    
    is_error <- FALSE
    tryCatch({
      ur<-uniroot(dummy,interval = c(0,up))$root
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

#' d>2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using a KDE approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
radial.thresh.KDE = function(r,w,tau=0.95,bww,bwr=0.05,up=30){ #ker.pdf=mv.Gaussian.ker.pdf, ker.cdf=Gaussian.ker.cdf){
  require(mvtnorm)
  
  # r, w      -> vector and matrix, angles defined by the L1 norm
  # bwr  -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)
  num.cols = dim(w)[2]
  
  r0w = apply(w, 1, function(ww){
    
    weightsw.i = lapply(1:ncol(bww),function(bww.idx){
      dnorm(x = w[,bww.idx],mean=ww[bww.idx], sd = bww[,bww.idx])
    })
    weightsw = apply(do.call(cbind,weightsw.i),1,prod)
    
    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
      # mean(weightsw*mean(rc<r))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,up))
    return(ur$root)
  })
  return(r0w)
}

#' d>2 setting: Evaluate a high threshold of R|W at a set of new angles W using a KDE method
#' 
#' @param wpts a vector of length n or a matrix with n rows and 2 columns
#' @param r Values of radial variable used to create the KDE threshold
#' @param w values of angular variable (on unit simplex) used to create the KDE threshold
#' @param tau level at which to calculate threshold
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
KDE.thresh.eval = function(wpts,r,w,tau=0.95,bww=NULL,bwr=0.05,up=30){
  require(mvtnorm)
  
  if(is.vector(wpts)){wpts = matrix(wpts,nrow=1)}
  
  n.dims = dim(w)[2]
  
  # if(is.null(bww)){
  #   bww = pdfCluster::kepdf(x = w, kernel = "gaussian", bwtype="adaptive",
  #                           alpha=0.5)@par$hx
  # } else if (length(bww==1)) {
  #   bww = matrix(bww,nrow(w),ncol(w))  # make the scalar value a matrix to be compatible with the following code
  # }
  
  r.tau.wpts = apply(wpts,1,function(wpts.i){
    if(sum(wpts.i)>1 | any(is.na(wpts.i)) | any(wpts.i<0)){
      return(NA)
    } else{
      
      weightsw.i = lapply(1:ncol(bww),function(bww.idx){
        dnorm(x = w[,bww.idx],mean=wpts.i[bww.idx], sd = bww[,bww.idx])
      })
      weightsw = apply(do.call(cbind,weightsw.i),1,prod)
      
      ccdf<-function(rc){
        mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
        # mean(weightsw*mean(rc<r))/mean(weightsw)
      }
      dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
      is_error <- FALSE
      tryCatch({
        ur<-uniroot(dummy,interval = c(0,up))$root
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

checkfn = function(quant,inp){
  inp*(quant - ifelse(inp<0,1,0))
}

###########################################################################


###############################################################################
################################  d=2  ########################################
###############################################################################

# Given angles and reference angles, return a list containing the angles and
# the indices of reference angles
which.adj.angles.2d = function(angles,locs,norm=NULL,marg="pos"){
  # w         -> a vector or 2-column matrix of angles in the d=1 unit simplex
  # locs      -> vector of length Nmesh, the reference angles in the d=1 unit simplex
  
  n.pars = length(locs)
  
  if(marg=="pos"){
    lst.val = lapply(angles,function(w){
      for(j in 1:(n.pars-1)){
        if(w>=locs[j] & w<=locs[j+1]){
          break
        }
      }
      return(list(w=w,
                  idx.locs=c(j,j+1)))
    })
  } else if(marg=="Rd"){
    lst.val = lapply(angles,function(w){
      for(j in 1:(n.pars)){
        idx.low = j
        idx.up = j+1
        val.low = locs[j]
        val.up = ifelse(idx.up>length(locs),-2,locs[j+1])
        if(w>=val.low & w<=val.up){
          break
        }
      }
      idx.up =  ifelse(idx.up>length(locs),1,idx.up)
      return(list(w=w,
                  idx.locs=c(idx.low,idx.up)))
    })
  }
  return(lst.val)
}

G.vol.2d = function(gauge.pars,par.locs,marg="pos"){
  if(marg=="pos"){
    val = 0.5*sum(sapply(c(1:(length(par.locs)-1)),
                         function(j){
                           par.locs.1.x = par.locs[j]
                           par.locs.2.x = par.locs[j+1]
                           par.locs.1.y = 1-par.locs[j]
                           par.locs.2.y = 1-par.locs[j+1]
                           return(gauge.pars[j]*gauge.pars[j+1]*
                                    abs(par.locs.1.x*par.locs.2.y - par.locs.2.x*par.locs.1.y))
                         }))
    return(val)
  } else if(marg=="Rd"){
    val = 0.5*sum(sapply(c(1:(length(par.locs))),
                         function(j){
                           loc1 = par.locs[j]
                           loc2 = ifelse(j+1>length(par.locs),-2,par.locs[j+1])
                           
                           g.pars.1 = gauge.pars[j]
                           g.pars.2 = ifelse(j+1>length(par.locs),gauge.pars[1],gauge.pars[j+1])
                           
                           par.locs.1 = pol2cart.L1Rd(w=loc1)
                           par.locs.2 = pol2cart.L1Rd(w=loc2)
                           par.locs.1.x = par.locs.1[1,1]
                           par.locs.2.x = par.locs.2[1,1]
                           par.locs.1.y = par.locs.1[1,2]
                           par.locs.2.y = par.locs.2[1,2]
                           return(g.pars.1*g.pars.2*
                                    abs(par.locs.1.x*par.locs.2.y - par.locs.2.x*par.locs.1.y))
                         }))
    return(val)
  }
}

###############################################################################
################################  d>2  ########################################
###############################################################################

which.adj.angles = function(angles,locs){
  # For each angle, where does it live in the partition of the simplex?
  
  require(geometry)
  
  # angles are d-dimensional
  
  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  if(is.null(dim(angles))){
    angles=t(as.matrix(angles))
  }
  
  num.cols = dim(locs)[2]
  
  del.tri = geometry::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  tsearchn.output = tsearchn(x=locs[,-num.cols],#rbind(locs,0)[,-num.cols],
                             t=del.tri$tri,
                             xi=matrix(angles[,-num.cols],ncol=num.cols-1))
  locs.idx = tsearchn.output$idx
  w.adj.angles = lapply(c(1:nrow(angles)),function(i){
    return(list(w=angles[i,],                                # the angle of interest, w
                tri.idx=locs.idx[i],                         # the index of the triangle where w belongs
                loc.idx=del.tri$tri[locs.idx[i],],           # the indices of the enclosing vertices
                vertices=locs[del.tri$tri[locs.idx[i],],]))  # the enclosing vertices (a matrix of size d x d)
  })
  return(w.adj.angles)
}

adj.DT = function(par.locs){
  # given reference angles, return the indices of the neighboring
  # triangles in the Delaunay triangulation
  
  num.cols = dim(par.locs)[2]
  del.tri = geometry::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)
  tri = del.tri$tri
  num.tris = nrow(tri)
  adj.lst = list(list(NULL))
  for(k in 1:num.tris){
    tri.k = tri[k,]
    adj.idx = NULL
    for(j in 1:num.tris){
      if(sum(table(c(tri.k,tri[j,]))==2)==(num.cols-1)){
        # if this is satisfied, tirangles i and j are neighbours
        adj.idx = c(adj.idx,j)
      }
    }
    adj.lst[[k]] = sort(adj.idx)
  }
  return(adj.lst)
}

ij.couples = function(par.locs){
  # given N par locs (nodes), for each node, find the (i,j) pairings in the penalty
  require(dplyr)
  
  num.cols = dim(par.locs)[2]
  tri = geometry::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)$tri
  adj.list = adj.DT(par.locs)
  
  ij.couples.lst = list(list(NULL))  # eventually will be of length N
  for(l in 1:nrow(par.locs)){
    #i=20
    # which triangle has location "i" as a vertex?
    which.has.l = apply(tri,1,function(ll){
      l %in% ll
    })
    if(sum(which.has.l)==1){
      # only 1 face, can't compare gradient to another
      ij.couples.lst[[l]] = NA
      next
    } else if(sum(which.has.l)>1){
      # which.grad.gs = (1:length(grad.g))[which.has.i]
      idx.w.l = c(1:length(adj.list))[which.has.l]
      adj.list.l = adj.list[idx.w.l]
      names(adj.list.l) = idx.w.l
      adj.list.l = lapply(adj.list.l, function(vec){
        return(vec[vec %in% idx.w.l])
      })
      ij.couples = lapply(1:length(idx.w.l), function(idx){
        idx.w.l.idx = idx.w.l[idx]
        lapply(adj.list.l[[idx]], function(vec){
          lapply(vec,function(vec.i){
            return(sort(c(idx.w.l.idx,vec.i)))
          })
        })
      })
      ij.couples = matrix(unlist(ij.couples),ncol=2,byrow=T)
      ij.couples = unname(as.matrix(distinct(as.data.frame(ij.couples,col.names=NULL))))
      ij.couples.lst[[l]] = ij.couples
    }
  }
  return(ij.couples.lst)
}

G.vol = function(gauge.pars,par.locs){
  num.cols = dim(par.locs)[2]
  
  # Method 2: coordinate geometry
  del.tri = geometry::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)
  nodes = del.tri$tri
  vol = (1/factorial(num.cols))*sum(apply(nodes,1,function(loc){
    pp = gauge.pars[loc]
    ll = par.locs[loc,]
    return(abs(det(ll*pp)))
  }))
  return(vol)
}

# Given angles and reference angles, return a list containing the reference angles,
# their indices, and how many angles are in each region
n.DT.region = function(which.adj.angles.res,locs){
  # which.adj.angles.res  -> output of which.adj.angles()
  # locs                  -> the reference angles in the d=2 unit simplex
  
  all.regions = lapply(which.adj.angles.res, function(lst) sort(lst$loc.idx))
  idx.locs = unique(all.regions)
  
  res = lapply(idx.locs,function(loc){
    freq = sum(sapply(which.adj.angles.res, function(lst) all(sort(lst$loc.idx) == loc)))
    return(list(reg = locs[loc,],  # the region is defined by these vertices
                loc = loc,
                freq=freq))        # there are this many angles in the region
  })
  
  return(res)
}

# # Given d-collumn exceedance angles, get reference angles such that each has
# # sufficient data to estmate
# get.par.logs = function(data){
#
#   res
#
#   return(res)
# }

###################################################################

# other functions

# sampling from f=e^-g

f.mcmc.g.2d<-function(niter,nburn,theta,alpha,thin,g){
  # initialize
  w<-rbeta(1,alpha,alpha)
  x<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)
  draws<-matrix(ncol=2,nrow=niter)
  it<-1
  while(it<=niter){
    w<-rbeta(1,alpha,alpha)
    xcan<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)
    
    accn<-g(xy=x,par=theta)*dbeta(x[1]/(x[1]+x[2]),alpha,alpha)*(xcan[1]+xcan[2])^2
    accd<-g(xy=xcan,par=theta)*dbeta(xcan[1]/(xcan[1]+xcan[2]),alpha,alpha)*(x[1]+x[2])^2
    
    if(runif(1)<accn/accd){x<-xcan}
    
    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

#########################################################


f.mcmc.g.3d<-function(niter,nburn,alpha=rep(1,3),thin,g){
  # initialize
  w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
  x<-rexp(1)*w/g(w)
  draws<-matrix(ncol=3,nrow=niter)
  it<-1
  while(it<=niter){
    w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
    xcan<-rexp(1)*w/g(w)
    
    accn<-g(x)*LaplacesDemon::ddirichlet(x/sum(x),alpha)*(sum(xcan))^3
    accd<-g(xcan)*LaplacesDemon::ddirichlet(xcan/sum(xcan),alpha)*(sum(x))^3
    
    if(runif(1)<accn/accd){x<-xcan}
    
    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

###########################################################

get_ref_angles = function(w, min.num = 30, gauss.corr=F){
  
  dd = dim(w)[2]
  
  # define par.locs
  # load("parlocs.RData")
  par.locs = expand.grid(replicate(dd-1, seq(0,1,by=0.2), simplify = FALSE))
  par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
  par.locs[,dd] = round(par.locs[,dd],3)
  par.locs = par.locs[par.locs[,dd]>=0,]
  par.locs = data.frame(rbind(diag(dd),as.matrix(unname(par.locs))))
  par.locs = as.matrix(unname(par.locs[!duplicated(par.locs),]))
  
  numm = rep(0,nrow(par.locs))
  while(any(numm < min.num)){
    
    nlocs = dim(par.locs)[1]
    which.adj.angles.res = which.adj.angles(w,par.locs)
    numm = sapply(1:nlocs, function(loc.fix){
      num = sum(sapply(which.adj.angles.res,function(lst){
        locs = lst$loc.idx
        return(loc.fix %in% locs)
      }))
      return(num)
    })
    # print(cbind(par.locs,numm))
    # print(numm)
    # which.rm = numm < min.num & c(1:nrow(par.locs)) > dd
    unique.min = min(unique(numm))
    which.rm = (numm < min.num) & (c(1:nrow(par.locs)) > dd) & (numm == unique.min)
    
    par.locs.new = par.locs[!which.rm,]
    if(nrow(par.locs.new) == nrow(par.locs)){
      ord = order(numm)
      which.rm = ord[ord>dd][1]
      par.locs = par.locs.new[-which.rm,]
    } else {
      par.locs = par.locs.new
    }
    
  }
  
  if(gauss.corr){
    par.locs[-c(1:dd),] = par.locs[-c(1:dd),] + rnorm(n=prod(dim(par.locs[-c(1:dd),])),sd=0.001)
    par.locs = ifelse(par.locs<0,0,par.locs)
    par.locs[-c(1:dd),] = par.locs[-c(1:dd),] / apply(par.locs[-c(1:dd),],1,sum)
  }
  
  return(par.locs)
}

get_ref_angles_5d = function(w,min.num=30,gauss.corr=F){
  
  # testing version - ignore for now.
  
  dd = dim(w)[2]
  
  # define par.locs
  # load("parlocs.RData")
  par.locs = expand.grid(replicate(dd-1, seq(0,1,by=0.1), simplify = FALSE))
  par.locs = cbind(par.locs,1-apply(par.locs,1,sum))
  par.locs[,dd] = round(par.locs[,dd],3)
  par.locs = par.locs[par.locs[,dd]>=0,]
  par.locs = data.frame(rbind(diag(dd),as.matrix(unname(par.locs))))
  par.locs = as.matrix(unname(par.locs[!duplicated(par.locs),]))
  
  # par.locs2 = rbind(diag(5),
  #                  rep(1/5,5),
  #                  c(1/4,1/4,1/4,1/4,0),c(1/4,1/4,1/4,0,1/4),c(1/4,1/4,0,1/4,1/4),c(1/4,0,1/4,1/4,1/4),c(0,1/4,1/4,1/4,1/4),
  #                  c(1/3,1/3,1/3,0,0),c(1/3,1/3,0,1/3,0),c(1/3,0,1/3,1/3,0),c(0,1/3,1/3,1/3,0),
  #                  c(1/3,1/3,0,0,1/3),c(1/3,0,1/3,0,1/3),c(0,1/3,1/3,0,1/3),
  #                  c(1/3,0,0,1/3,1/3), c(0,1/3,0,1/3,1/3), c(0,0,1/3,1/3,1/3),
  #                  c(0.5,0.5,0,0,0),c(0.5,0,0.5,0,0),c(0.5,0,0,0.5,0),c(0.5,0,0,0,0.5),
  #                  c(0,0.5,0.5,0,0),c(0,0.5,0,0.5,0),c(0,0.5,0,0,0.5),
  #                  c(0,0,0.5,0.5,0),c(0,0,0.5,0,0.5),c(0,0,0,0.5,0.5)
  # )
  # par.locs = rbind(diag(dd),par.locs)
  # par.locs = as.matrix(unname(par.locs[!duplicated(data.frame(par.locs)),]))
  
  numm = rep(0,nrow(par.locs))
  while(any(numm < min.num)){
    
    nlocs = dim(par.locs)[1]
    which.adj.angles.res = which.adj.angles(w,par.locs)
    numm = sapply(1:nlocs, function(loc.fix){
      num = sum(sapply(which.adj.angles.res,function(lst){
        locs = lst$loc.idx
        return(loc.fix %in% locs)
      }))
      return(num)
    })
    # print(cbind(par.locs,numm))
    # print(numm)
    unique.min = min(unique(numm))
    which.rm = (numm < min.num) & (c(1:nrow(par.locs)) > dd) & (numm == unique.min)
    par.locs.new = par.locs[!which.rm,]
    if(nrow(par.locs.new) == nrow(par.locs)){
      ord = order(numm)
      which.rm = ord[ord>dd][1]
      par.locs = par.locs.new[-which.rm,]
    } else {
      par.locs = par.locs.new
    }
    
  }
  
  if(gauss.corr){
    par.locs[-c(1:dd),] = par.locs[-c(1:dd),] + rnorm(n=prod(dim(par.locs[-c(1:dd),])),sd=0.001)
    par.locs = ifelse(par.locs<0,0,par.locs)
    par.locs[-c(1:dd),] = par.locs[-c(1:dd),] / apply(par.locs[-c(1:dd),],1,sum)
  }
  
  return(par.locs)
}



# d=2 pwl gauge function def

# Given angles and reference angles, return a list containing the angles and
# the indices of reference angles
which.adj.angles.2d = function(angles,locs,norm=NULL,marg="pos"){
  # w         -> a vector or 2-column matrix of angles in the d=1 unit simplex
  # locs      -> vector of length Nmesh, the reference angles in the d=1 unit simplex

  n.pars = length(locs)

  if(marg=="pos"){
    lst.val = lapply(angles,function(w){
      for(j in 1:(n.pars-1)){
        if(w>=locs[j] & w<=locs[j+1]){
          break
        }
      }
      return(list(w=w,
                  idx.locs=c(j,j+1)))
    })
  } else if(marg=="Rd"){
    lst.val = lapply(angles,function(w){
      for(j in 1:(n.pars)){
        idx.low = j
        idx.up = j+1
        val.low = locs[j]
        val.up = ifelse(idx.up>length(locs),-2,locs[j+1])
        if(w>=val.low & w<=val.up){
          break
        }
      }
      idx.up =  ifelse(idx.up>length(locs),1,idx.up)
      return(list(w=w,
                  idx.locs=c(idx.low,idx.up)))
    })
  }
  return(lst.val)
}


# given 3 parameters (par) located at 3 reference angles (par.locs),
# return the gauge function value at the point xyz
gfun.simple.d2pwlin.L1 = function(w,par,par.locs){
  # xyz      =
  # par      = vector of length 2
  # par.locs = vector of length 2

  r1 = par[1]
  r2 = par[2]

  w.sort = sort(par.locs)
  w.low = w.sort[1]
  w.up = w.sort[2]

  dx = r2*w.up - r1*w.low
  dy = r2*(1-w.up) - r1*(1-w.low)
  a = dy/dx
  b = r1*(1-w.low - a*w.low)
  if(dx==0){
    return(w/(r1*w.low))
  }
  else{
    return((1-w-(a*w))/b)
  }
}


pwlin.g.vals.2d = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.

  sapply(w.adj.angles, function(lst){
    which.angles = lst$idx.locs
    return(gfun.simple.d2pwlin.L1(w=lst$w,par=par[which.angles],par.locs=par.locs[which.angles]))
  })
}

gfun.2d = function(x, par, ref.angles){
  # x-> matrix of 2 columns
  if(is.null(dim(x)) & length(x)==1){
    x = cbind(x,1-x)
  }
  if(is.null(dim(x)) & length(x)>1){
    x = matrix(x,nrow=1)
  }

  rad.inp = apply(x,1,sum)
  angle.inp = x[,1] / rad.inp

  w.adj.angles = which.adj.angles.2d(angles=angle.inp, locs=ref.angles)

  gval = rad.inp*pwlin.g.vals.2d(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)

  return(gval)
}

######################################################################


# plot the 3-d projection of the gauge function

gfun.simple.pwlin = function(xyz,par,par.locs){
  # par      = in R^3_+
  # par.locs = angles (in the simplex) where the parameters are location

  num.cols = dim(par.locs)[2]

  coplanar.mat = do.call(rbind,lapply(c(1:num.cols)[-1],function(i){
    (par[1]*par.locs[1,])-(par[i]*par.locs[i,])
  }))
  norm.vec = suppressWarnings({
    c(1,-1)*sapply(c(1:num.cols),function(i){det(coplanar.mat[,-i])})
  })
  sum(norm.vec * xyz) / sum(norm.vec * par[1] * par.locs[1,])
}

# evaluating the gauge at a set of angles, returns a vector
pwlin.g.vals = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.

  sapply(w.adj.angles, function(lst){
    which.angles = lst$loc.idx
    if(any(is.na(lst$w)) | all(which.angles==0) | any(lst$w<0) | any(is.na(which.angles))){  # if angle is not in the positive orthant OR if angle lies on line
      return(NA)
    } else {
      return(gfun.simple.pwlin(xyz=lst$w,par=par[which.angles],par.locs=par.locs[which.angles,]))
    }
  })
}

# qhull.options <- function(options, output.options, supported_output.options, full=FALSE) {
#   if (full) {
#     if (!is.null(output.options)) {
#       stop("full and output.options should not be specified together")
#     }
#     output.options = TRUE
#     ## Enable message in 0.4.1
#     ## Turn to warning in 0.4.7
#     message("delaunayn: \"full\" option is deprecated; adding \"Fa\" and \"Fn\" to options.
#       Please update your code to use \"output.options=TRUE\" or set \"output.options\" to a
#       string containing desired QHull options.")
#   }
#   
#   if (is.null(output.options)) {
#     output.options <- ""
#   }
#   if (is.logical(output.options)) {
#     if (output.options) {
#       output.options <- paste(supported_output.options, collapse=" ")
#     } else {
#       output.options  <- ""
#     }
#   }
#   if (!is.character(output.options)) {
#     stop("output.options must be a string, logical or NULL")
#   }
#   
#   ## Input sanitisation
#   options <- paste(options, output.options, collapse=" ")
#   return(options)
# }
# 
# delaunayn <-
#   function(p, options=NULL, output.options=NULL, full=FALSE) {
#     tmp_stdout <- tempfile("Rf")
#     tmp_stderr <- tempfile("Rf")
#     on.exit(unlink(c(tmp_stdout, tmp_stderr)))
#     
#     ## Coerce the input to be matrix
#     if (is.data.frame(p)) {
#       p <- as.matrix(p)
#     }
#     
#     ## Make sure we have real-valued input
#     storage.mode(p) <- "double"
#     
#     ## We need to check for NAs in the input, as these will crash the C
#     ## code.
#     if (any(is.na(p))) {
#       stop("The first argument should not contain any NAs")
#     }
#     
#     ## Default options
#     default.options <- "Qt Qc Qx"
#     if (ncol(p) < 4) {
#       default.options <- "Qt Qc Qz"
#     }
#     if (is.null(options)) {
#       options <- default.options
#     }
#     
#     ## Combine and check options
#     options <- tryCatch(qhull.options(options, output.options, supported_output.options  <- c("Fa", "Fn"), full=full), error=function(e) {stop(e)})
#     
#     ## It is essential that delaunayn is called with either the QJ or Qt
#     ## option. Otherwise it may return a non-triangulated structure, i.e
#     ## one with more than dim+1 points per structure, where dim is the
#     ## dimension in which the points p reside.
#     if (!grepl("Qt", options) & !grepl("QJ", options)) {
#       options <- paste(options, "Qt")
#     }
#     
#     out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout, tmp_stderr, PACKAGE="geometry")
#     
#     ## Check for points missing from triangulation, but not in the case
#     ## of a degenerate trianguation (zero rows in output)
#     if (nrow(out$tri) > 0) {
#       missing.points <- length(setdiff(seq(1,nrow(p)), unique(as.vector(out$tri))))
#       if (missing.points > 0) {
#         warning(paste0(missing.points, " points missing from triangulation.
# It is possible that setting the 'options' argument of delaunayn may help.
# For example:
# options = \"", default.options, " Qbb\"
# options = \"", default.options, " QbB\"
# If these options do not work, try shifting the centre of the points
# to the origin by subtracting the mean coordinates from every point."))
#       }
#     }
#     
#     # Remove NULL elements
#     out[which(sapply(out, is.null))] <- NULL
#     if (is.null(out$areas) & is.null(out$neighbours)) {
#       attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
#       return(out$tri)
#     }
#     class(out) <- "delaunayn"
#     out$p <- p
#     return(out)
#   }

# 
# which.adj.angles = function(angles,locs){
#   # For each angle, where does it live in the partition of the simplex?
# 
#   require(geometry)
# 
#   # angles are d-dimensional
# 
#   if("data.frame" %in% class(locs)){
#     locs=as.matrix(locs)
#   }
#   if("data.frame" %in% class(angles)){
#     angles=as.matrix(angles)
#   }
#   if(is.null(dim(angles))){
#     angles=t(as.matrix(angles))
#   }
# 
#   num.cols = dim(locs)[2]
# 
#   del.tri = delaunayn(p=locs[,-num.cols], output.options=TRUE)
#   tsearchn.output = tsearchn(x=locs[,-num.cols],#rbind(locs,0)[,-num.cols],
#                              t=del.tri$tri,
#                              xi=matrix(angles[,-num.cols],ncol=num.cols-1))
#   locs.idx = tsearchn.output$idx
#   w.adj.angles = lapply(c(1:nrow(angles)),function(i){
#     return(list(w=angles[i,],                                # the angle of interest, w
#                 tri.idx=locs.idx[i],                         # the index of the triangle where w belongs
#                 loc.idx=del.tri$tri[locs.idx[i],],           # the indices of the enclosing vertices
#                 vertices=locs[del.tri$tri[locs.idx[i],],]))  # the enclosing vertices (a matrix of size d x d)
#   })
#   return(w.adj.angles)
# }
gfun.pwl = function(x, par, ref.angles){
  # x-> matrix of d columns
  if(is.null(dim(x))){
    x = matrix(x,nrow=1)
  }

  num.cols = dim(x)[2]

  rad.inp = apply(x,1,sum)
  angle.inp = x / rad.inp
  is.valid = apply(angle.inp,1,function(vec) all(!is.na(vec)))
  w.adj.angles = replicate(nrow(angle.inp),
                           list(w=rep(NA,num.cols),
                                loc.idx=rep(NA,num.cols),
                                vertices=matrix(NA,num.cols,num.cols)), FALSE)
  inn = as.matrix(angle.inp[is.valid,],ncol=num.cols,byrow=T)
  if(ncol(inn)==1){
    inn = t(inn)
  }
  w.adj.angles[is.valid] = which.adj.angles(angles=inn, locs=ref.angles)

  gval = rad.inp*pwlin.g.vals(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)
  return(gval)
}

proj.g.fn = function(gfun,w,which.w,nm,...){
  # gfun -> gauge that takes in 4-dim vectors
  # w -> 3-min input
  # which.w -> which index to take min over

  w.inp = matrix(NA,nrow=nm,ncol=4)
  w.inp[,-which.w] = matrix(as.numeric(w),ncol=3,nrow=nm,byrow=T)
  w.inp[,which.w] = seq(0,1.3,length.out=nm)

  return(min(gfun(w.inp,...)))
}

f.mcmc.g.2d<-function(niter,nburn,theta,alpha,thin,g){
  # initialize
  w<-rbeta(1,alpha,alpha)
  x<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)
  draws<-matrix(ncol=2,nrow=niter)
  it<-1
  while(it<=niter){
    w<-rbeta(1,alpha,alpha)
    xcan<-rexp(1)*c(w,1-w)/g(c(w,1-w),par=theta)
    
    accn<-g(xy=x,par=theta)*dbeta(x[1]/(x[1]+x[2]),alpha,alpha)*(xcan[1]+xcan[2])^2
    accd<-g(xy=xcan,par=theta)*dbeta(xcan[1]/(xcan[1]+xcan[2]),alpha,alpha)*(x[1]+x[2])^2
    
    if(runif(1)<accn/accd){x<-xcan}
    
    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

#########################################################


f.mcmc.g.3d<-function(niter,nburn,alpha=rep(1,3),thin,g){
  # initialize
  w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
  x<-rexp(1)*w/g(w)
  draws<-matrix(ncol=3,nrow=niter)
  it<-1
  while(it<=niter){
    w<-as.numeric(LaplacesDemon::rdirichlet(1,alpha))
    xcan<-rexp(1)*w/g(w)
    
    accn<-g(x)*LaplacesDemon::ddirichlet(x/sum(x),alpha)*(sum(xcan))^3
    accd<-g(xcan)*LaplacesDemon::ddirichlet(xcan/sum(xcan),alpha)*(sum(x))^3
    
    if(runif(1)<accn/accd){x<-xcan}
    
    draws[it,]<-x
    it<-it+1
  }
  return(draws[-(1:nburn),][seq(1,(niter-nburn),by=thin),])
}

############################################################################

# Extra, manuscript-specific sampling functions

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
    iw <- iweights.pwl(k = k, r0w = r0w, w = w, gfun = gfun, par = par, shape=shape)
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

