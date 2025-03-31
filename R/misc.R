# A clone of the geometry::delaunayn function, without the recent
# update (Feb 2025) with typos
# delaunayn = function (p, options = NULL, output.options = NULL, full = FALSE) {
#   tmp_stdout <- tempfile("Rf")
#   tmp_stderr <- tempfile("Rf")
#   on.exit(unlink(c(tmp_stdout, tmp_stderr)))
#   if (is.data.frame(p)) {
#     p <- as.matrix(p)
#   }
#   storage.mode(p) <- "double"
#   if (any(is.na(p))) {
#     stop("The first argument should not contain any NAs")
#   }
#   if (is.null(options)) {
#     if (ncol(p) < 4) {
#       options <- "Qt Qc Qz"
#     }
#     else {
#       options <- "Qt Qc Qx"
#     }
#   }
#   options <- tryCatch(qhull.options(options, output.options,
#                                     supported_output.options <- c("Fa", "Fn"), full = full),
#                       error = function(e) {
#                         stop(e)
#                       })
#   if (!grepl("Qt", options) & !grepl("QJ", options)) {
#     options <- paste(options, "Qt")
#   }
#   out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout,
#                tmp_stderr, PACKAGE = "geometry")
#   out[which(sapply(out, is.null))] <- NULL
#   if (is.null(out$areas) & is.null(out$neighbours)) {
#     attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
#     return(out$tri)
#   }
#   class(out) <- "delaunayn"
#   out$p <- p
#   return(out)
# }

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
    options <- tryCatch(PWLExtremes::qhull.options(options, output.options, supported_output.options  <- c("Fa", "Fn"), full=full), error=function(e) {stop(e)})

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
  require(geometry)

  # angles are d-dimensional

  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  if(is.null(dim(angles))){
    angles=as.matrix(angles,nrow=1)
  }

  num.cols = dim(locs)[2]

  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  w.adj.angles = tsearchn(x=locs[,-num.cols],#rbind(locs,0)[,-num.cols],
                          t=del.tri$tri,
                          xi=matrix(angles[,-num.cols],ncol=num.cols-1))
  locs.idx =w.adj.angles$idx
  w.adj.angles = lapply(c(1:nrow(angles)),function(i){
    return(list(w=angles[i,],                                           # the angle
                loc.idx=del.tri$tri[locs.idx[i],],
                vertices=locs[del.tri$tri[locs.idx[i],],]))  # the enclosing vertices
  })
  return(w.adj.angles)
}

G.vol = function(gauge.pars,par.locs){
  num.cols = dim(par.locs)[2]

  # Method 2: coordinate geometry
  del.tri = PWLExtremes::delaunayn(p=par.locs[,-num.cols], output.options=TRUE)
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

# cart2pol.L1Rd = function(x){
#   if(dim(x)[2]!=2){
#     stop("input needs to be a 2 column matrix")
#   }
#   r<-apply(x,1,function(xx) sum(abs(xx)))
#   w.cart<-x/r
#   eps.val = ifelse(w.cart[,2]>=0,1,-1)
#   w.pol=(1-w.cart[,1])*sign(w.cart[,2])
#   # w.pol[r==0] = 0
#   return(list(r=r,w=w.pol))
# }
#
# cos1 = function(q){
#   q = q %% 4
#   q[q>2] = q[q>2]-4
#   return(1-abs(q))
# }
# sin1 = function(q){
#   return(cos1(q-1))
# }
# pol2cart.L1Rd = function(r=NULL,w){
#   w1 = cos1(w)
#   w2 = cos1(w-1)
#   if(is.null(r)){
#     return(cbind(w1,w2))
#   } else {
#     return(r*cbind(w1,w2))
#   }
# }

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

