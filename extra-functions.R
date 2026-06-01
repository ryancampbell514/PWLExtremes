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
      weightsw = dmvnorm(w[,-n.dims],mean=wpts.i[-n.dims],sigma=(bww^2)*diag(n.dims-1))
      
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
                                      bww=thresh.fit$bww)
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
    w = w[,-ncol(w)]
    require(pdfCluster)
    if(is.null(bww)){
      require(pdfCluster)
      bww = as.numeric(pdfCluster::kepdf(x = w, kernel = "gaussian", bwtype="adaptive",
                                         alpha=alpha)@par$hx)  # can play with the apha hyperparameter
    }
    r0w = radial.thresh.KDE.2d(r=r,w=w,tau=tau,bww=bww,up=up)
  } else if(dim(w)[2]>2 && method=="KDE") {
    w = w[,-ncol(w)]
    if(is.null(bww)){
      require(pdfCluster)
      bww = pdfCluster::kepdf(x = w, kernel = "gaussian", bwtype="adaptive",
                              alpha=alpha)@par$hx
    } else if (length(bww)==1) {
      bww = matrix(bww,nrow(w),ncol(w))  # make the scalar value a matrix to be compatible with the code in radial.thresh.KDE
    } else {stop("bww must be NULL or a single numeric value")}
    r0w = radial.thresh.KDE(r=r,w=w,tau=tau,bww=bww,up=up)
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

# Empirical threshold estimation

#' Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  stop("General d-dimension empirical threshold estimation not yet implemented. Check back later.")
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  qu.r <- array(0, dim = rep(bin.mesh, (dim(w)[2]-1)))
  array.idx = expand.grid(replicate(dim(w)[2]-1, c(1:bin.mesh), simplify = FALSE))
  for(rrow in 1:nrow(array.idx)){
    # TODO: SPEED UP THE NEXT 4 LINES OF CODE
    idx = array.idx[rrow,]
    ind = rep(TRUE,nrow(w))
    for(ii in 1:length(idx)){
      ind = ind & (wsl[as.numeric(idx[ii])] < w[,ii]) & (w[,ii] < wsu[as.numeric(idx[ii])])
    }
    r1 <- r[ind]
    
    if (sum(ind) > 2) {
      qu.r[matrix(as.numeric(idx),1)] <- quantile(r1, tau)
    }
    else {
      qu.r[matrix(as.numeric(idx),1)] <- NA
    }
  }
  
  r0w <- rep(0,nrow(w))
  # array.idx = expand.grid(replicate(dim(w)[2]-1, c(1:bin.mesh), simplify = FALSE))
  for(i in 1:length(r0w)){
    indArray <- array(NA, dim = rep(bin.mesh,dim(w)[2]-1))
    for(rrow in 1:nrow(array.idx)){
      # TODO: SPEED UP THE NEXT 4 LINES OF CODE
      idx = array.idx[rrow,]
      ind = TRUE
      for(ii in 1:length(idx)){
        ind = ind & (wsl[as.numeric(idx[ii])] < w[i,ii]) & (w[i,ii] < wsu[as.numeric(idx[ii])])
      }
      indArray[matrix(as.numeric(idx),1)] = ind
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  return(r0w)
}

#' d=2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.2d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(dim(w))){
    stop("w needs to be a n x 2 matrix")
  }
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w = w[,1]  # make into a vector
  
  qu.r <- NULL
  for (j in 1:bin.mesh) {
    re <- r[wsl[j] < w & w < wsu[j]]
    qu.r[j] <- quantile(re, tau)
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indvec <- rep(NA, bin.mesh)
    for (j in 1:bin.mesh) {
      indvec[j] <- wsl[j] < w[i] & w[i] < wsu[j]
    }
    r0w[i] <- mean(qu.r[indvec], na.rm = T)
  }
  return(list(r0w=r0w,             # the threshold
              quant.grid=qu.r))
}

#' d=2 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.2d = function(wpts,quant.grid,bin.mesh,overlap){
  if(is.null(dim(wpts))){
    stop("wpts needs to be a n x 2 matrix")
  }
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w = wpts[,1]  # make into a vector
  
  n <- length(w)
  r0w <- numeric(n)
  for (i in 1:n) {
    indvec <- rep(NA, bin.mesh)
    for (j in 1:bin.mesh) {
      indvec[j] <- wsl[j] < w[i] & w[i] < wsu[j]
    }
    r0w[i] <- mean(quant.grid[indvec], na.rm = T)
  }
  
  return(r0w)
}

#' d=3 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.3d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[,1]  # make into a vector
  w2 = w[,2]  # make into a vector
  
  qu.r <- matrix(0, bin.mesh, bin.mesh)
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      ind <- 
        wsl[j] < w1 & w1 < wsu[j] & 
        wsl[k] < w2 & w2 < wsu[k]
      r1 <- r[ind]
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
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indMat <- matrix(NA, bin.mesh, bin.mesh)
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        indMat[j, k] <- 
          wsl[j] < w1[i] & w1[i] < wsu[j] &
          wsl[k] < w2[i] & w2[i] < wsu[k]
      }
    }
    r0w[i] <- mean(qu.r[indMat], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=3 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.3d = function(wpts,quant.grid,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[,1]  # make into a vector
  w2 = wpts[,2]  # make into a vector
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indMat <- matrix(NA, bin.mesh, bin.mesh)
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        indMat[j, k] <- 
          wsl[j] < w1[i] & w1[i] < wsu[j] &
          wsl[k] < w2[i] & w2[i] < wsu[k]
      }
    }
    r0w[i] <- mean(quant.grid[indMat], na.rm = T)
  }
  
  return(r0w)
}

#' d=4 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.4d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[, 1]
  w2 = w[, 2]
  w3 = w[, 3]
  
  qu.r <- array(0, dim = rep(bin.mesh,3))
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      for (l in 1:bin.mesh) {
        ind <- 
          wsl[j] < w1 & w1 < wsu[j] & 
          wsl[k] < w2 & w2 < wsu[k] & 
          wsl[l] < w3 & w3 < wsu[l]
        r1 <- r[ind]
        if (sum(ind) > 10) {
          qu.r[j, k, l] <- quantile(r1, tau)
        }
        else {
          qu.r[j, k, l] <- NA
        }
      }
    }
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,3))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          indArray[j, k, l] <- 
            wsl[j] < w1[i] & w1[i] < wsu[j] & 
            wsl[k] < w2[i] & w2[i] < wsu[k] &
            wsl[l] < w3[i] & w3[i] < wsu[l]
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=4 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.4d = function(wpts,quant.grid, bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[, 1]
  w2 = wpts[, 2]
  w3 = wpts[, 3]
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,3))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          indArray[j, k, l] <- 
            wsl[j] < w1[i] & w1[i] < wsu[j] & 
            wsl[k] < w2[i] & w2[i] < wsu[k] &
            wsl[l] < w3[i] & w3[i] < wsu[l]
        }
      }
    }
    r0w[i] <- mean(quant.grid[indArray], na.rm = T)
  }
  
  return(r0w)
}

#' d=5 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.5d = function(r, w, tau = 0.95,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[, 1]
  w2 = w[, 2]
  w3 = w[, 3]
  w4 = w[, 4]
  
  qu.r <- array(0, dim = rep(bin.mesh,4))
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      for (l in 1:bin.mesh) {
        for(m in 1:bin.mesh) {
          ind <- 
            (wsl[j] < w1) & (w1 < wsu[j]) &
            (wsl[k] < w2) &( w2 < wsu[k]) &
            (wsl[l] < w3) &( w3 < wsu[l]) &
            (wsl[m] < w4) & (w4 < wsu[m])
          # print(c(,sum(ind)))
          # print(ind)
          r1 <- r[ind]
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
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,4))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          for (m in 1:bin.mesh) {
            indArray[j, k, l, m] <-
              wsl[j] < w1[i] & w1[i] < wsu[j] &
              wsl[k] < w2[i] & w2[i] < wsu[k] &
              wsl[l] < w3[i] & w3[i] < wsu[l] &
              wsl[m] < w4[i] & w4[i] < wsu[m]
          }
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=5 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.5d = function(wpts,quant.grid,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[, 1]
  w2 = wpts[, 2]
  w3 = wpts[, 3]
  w4 = wpts[, 4]
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,4))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          for (m in 1:bin.mesh) {
            indArray[j, k, l, m] <-
              wsl[j] < w1[i] & w1[i] < wsu[j] &
              wsl[k] < w2[i] & w2[i] < wsu[k] &
              wsl[l] < w3[i] & w3[i] < wsu[l] &
              wsl[m] < w4[i] & w4[i] < wsu[m]
          }
        }
      }
    }
    r0w[i] <- mean(quant.grid[indArray], na.rm = T)
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

qhull.options <- function(options, output.options, supported_output.options, full=FALSE) {
  if (full) {
    if (!is.null(output.options)) {
      stop("full and output.options should not be specified together")
    }
    output.options = TRUE
    ## Enable message in 0.4.1
    ## Turn to warning in 0.4.7
    message("delaunayn: \"full\" option is deprecated; adding \"Fa\" and \"Fn\" to options.
      Please update your code to use \"output.options=TRUE\" or set \"output.options\" to a
      string containing desired QHull options.")
  }
  
  if (is.null(output.options)) {
    output.options <- ""
  }
  if (is.logical(output.options)) {
    if (output.options) {
      output.options <- paste(supported_output.options, collapse=" ")
    } else {
      output.options  <- ""
    }
  }
  if (!is.character(output.options)) {
    stop("output.options must be a string, logical or NULL")
  }
  
  ## Input sanitisation
  options <- paste(options, output.options, collapse=" ")
  return(options)
}

delaunayn <-
  function(p, options=NULL, output.options=NULL, full=FALSE) {
    tmp_stdout <- tempfile("Rf")
    tmp_stderr <- tempfile("Rf")
    on.exit(unlink(c(tmp_stdout, tmp_stderr)))
    
    ## Coerce the input to be matrix
    if (is.data.frame(p)) {
      p <- as.matrix(p)
    }
    
    ## Make sure we have real-valued input
    storage.mode(p) <- "double"
    
    ## We need to check for NAs in the input, as these will crash the C
    ## code.
    if (any(is.na(p))) {
      stop("The first argument should not contain any NAs")
    }
    
    ## Default options
    default.options <- "Qt Qc Qx"
    if (ncol(p) < 4) {
      default.options <- "Qt Qc Qz"
    }
    if (is.null(options)) {
      options <- default.options
    }
    
    ## Combine and check options
    options <- tryCatch(qhull.options(options, output.options, supported_output.options  <- c("Fa", "Fn"), full=full), error=function(e) {stop(e)})
    
    ## It is essential that delaunayn is called with either the QJ or Qt
    ## option. Otherwise it may return a non-triangulated structure, i.e
    ## one with more than dim+1 points per structure, where dim is the
    ## dimension in which the points p reside.
    if (!grepl("Qt", options) & !grepl("QJ", options)) {
      options <- paste(options, "Qt")
    }
    
    out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout, tmp_stderr, PACKAGE="geometry")
    
    ## Check for points missing from triangulation, but not in the case
    ## of a degenerate trianguation (zero rows in output)
    if (nrow(out$tri) > 0) {
      missing.points <- length(setdiff(seq(1,nrow(p)), unique(as.vector(out$tri))))
      if (missing.points > 0) {
        warning(paste0(missing.points, " points missing from triangulation.
It is possible that setting the 'options' argument of delaunayn may help.
For example:
options = \"", default.options, " Qbb\"
options = \"", default.options, " QbB\"
If these options do not work, try shifting the centre of the points
to the origin by subtracting the mean coordinates from every point."))
      }
    }
    
    # Remove NULL elements
    out[which(sapply(out, is.null))] <- NULL
    if (is.null(out$areas) & is.null(out$neighbours)) {
      attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
      return(out$tri)
    }
    class(out) <- "delaunayn"
    out$p <- p
    return(out)
  }


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
  
  del.tri = delaunayn(p=locs[,-num.cols], output.options=TRUE)
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
