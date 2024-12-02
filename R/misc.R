
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




# # Angular density:
# # for now, unnormalised, need to numerically integrate over the mesh
# fW.2d = function(w,par,locs=NULL,true.gauge=NULL){
#   # w -> angle on the unit simplex
#   # par -> gauge params
#   # locs -> location of reference angles on the unit simplex
#   
#   if(is.null(locs) & is.null(true.gauge)){
#     stop("Specify reference angles or true gauge")
#   }
# 
#   # 1) put cartesian to simplex
#   if(length(w)==2 & sum(w)==1){
#     w = w[1]
#   }
#   
#   # 2) evaluate gauge at w
#   if(is.null(true.gauge)){
#     adj.angles = which.adj.angles.2d(angles=w, locs)
#     val = pwlin.g.vals.2d(adj.angles,par,locs)
#   } else {
#     val = true.gauge(c(w,1-w),par)
#   }
#   
#   
#   # 3) return value
#   return(val^(-2))
# }
# 
G.vol.2d = function(gauge.pars,par.locs,marg="pos"){
  if(marg=="pos"){
    val = 0.5*sum(sapply(c(1:(length(par.locs)-1)),
                   function(j){
                     # m = rbind(rep(0,2),
                     #           gauge.pars[j]*c(par.locs[j],1-par.locs[j]),
                     #           gauge.pars[j+1]*c(par.locs[j+1],1-par.locs[j+1]))
                     # m = cbind(m,1)
                     # return(det(m))
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
                       # loc2 = par.locs[j+1]
                       loc2 = ifelse(j+1>length(par.locs),-2,par.locs[j+1])
                       # print(loc2)
                       
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
################################  d=3  ########################################
###############################################################################

cor2cov <- function(V, sd) {
  V * tcrossprod(sd)
}


# Given angles and reference angles, return a list containing the angles and
# the indices of reference angles
which.adj.angles.3d = function(angles,locs){
  # w         -> a 2 or 3-column matrix of angles in the d=2 unit simplex
  # locs      -> the reference angles in the d=2 unit simplex
  
  if(is.null(dim(locs))){
    locs=matrix(locs,nrow=1)
  }
  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if(dim(locs)[2]==3){
      locs = locs[,-3]
  }
  # if(!("data.frame" %in% class(locs))){
  #   locs = data.frame(locs,row.names=NULL)
  #   names(locs) = c("x","y")
  # }
  # del.tri = interp::tri.mesh(locs) 
  del.tri = geometry::delaunayn(p=locs, output.options=TRUE)
  
  if(is.null(dim(angles))){
    angles=matrix(angles,nrow=1)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  if(dim(angles)[2]==3){
    angles=angles[,-3]
  }
  if(is.null(dim(angles))){
    angles=matrix(angles,nrow=1)
  }
  
  val = tsearchn(x=locs,t=del.tri$tri,xi=angles)
  # val$idx
  lst.val = lapply(1:length(val$idx),function(ii){
    return(list(w=c(angles[ii,],1-sum(angles[ii,])),
                idx.locs=del.tri$tri[val$idx[ii],]))
  })
  # lst.val = apply(angles,1,function(w){
  #   where.pt = interp::tri.find(tri.obj=del.tri,x=w[1],y=w[2])
  #   idx.locs = c(where.pt$i1,where.pt$i2,where.pt$i3)
  #   return(list(w=w,
  #               idx.locs=idx.locs))
  # })
  return(lst.val)
}

which.adj.angles.3d.lap = function(angles,locs){
  # w         -> a 2 or 3-column matrix of angles in the d=2 unit simplex
  # locs      -> the reference angles in the d=2 unit simplex
  if(dim(locs)[2]==3){
    # project to theta and phi
    locs = geometry::cart2sph(locs)[,-3]
  }
  if(!("data.frame" %in% class(locs))){
    locs = data.frame(locs,row.names=NULL)
    names(locs) = c("x","y")
  }
  stop("remove interp dependency, replace with geometry::delaunayn")
  del.tri = interp::tri.mesh(locs) 
  # plot.triSht(del.tri,do.labels=T)
  if(is.null(dim(angles))){
    angles=matrix(angles,nrow=1)
  }
  # print(angles)
  angles = geometry::cart2sph(angles)[,-3]
  # points(angles)
  lst.val = apply(angles,1,function(w){
    stop("remove interp dependency, replace with geometry::delaunayn")
    # where.pt = interp::tri.find(tri.obj=del.tri,x=w[1],y=w[2])
    idx.locs = c(where.pt$i1,where.pt$i2,where.pt$i3)
    return(list(w=as.numeric(geometry::sph2cart(matrix(c(w,1),nrow=1))),
                idx.locs=idx.locs))
  })
  return(lst.val)
}

# Given angles and reference angles, return a list containing the reference angles,
# their indices, and how many angles are in each region
n.DT.region.3d = function(which.adj.angles.res,locs){
  # which.adj.angles.res  -> output of which.adj.angles()
  # locs                  -> the reference angles in the d=2 unit simplex
  
  # if(dim(locs)[2]==3){
  #   locs = locs[,-3]
  # }
  # if(class(locs) != "data.frame"){
  #   locs = data.frame(locs,row.names=NULL)
  #   names(locs) = c("x","y")
  # }
  
  all.regions = lapply(which.adj.angles.res, function(lst) lst$idx.locs == loc)
  idx.locs = unique(all.regions)
  
  res = lapply(idx.locs,function(loc){
    freq = sum(sapply(which.adj.angles.res, function(lst) all(lst$idx.locs == loc)))
    return(list(reg = locs[loc,],
                freq=freq))
  })
}
# 
# # Angular density:
# # for now, unnormalised, need to numerically integrate over the mesh
# fW.3d = function(w,par,locs=NULL,true.gauge=NULL){
#   # w -> angle on the unit simplex
#   # par -> gauge params
#   # locs -> location of reference angles on the unit simplex
#   
#   if(is.null(locs) & is.null(true.gauge)){
#     stop("Specify reference angles or true gauge")
#   }
#   
#   # # 1) put to cartesian
#   if(length(w)==2){
#     w = cbind(w,1-sum(w))
#   }
#    
#   # 2) evaluate gauge at w
#   if(is.null(true.gauge)){
#     adj.angles = which.adj.angles.2d(angles=w, locs)
#     val = pwlin.g.vals.2d(adj.angles,par,locs)
#   } else {
#     val = true.gauge(w,par)
#   }
#   
#   # 3) return value
#   return(val^(-3))
# }

G.vol.3d = function(gauge.pars,par.locs){

  # par.locs.df = data.frame(par.locs[,-3],row.names=NULL)
  # names(par.locs.df) = c("x","y")
  # del.tri = interp::tri.mesh(par.locs.df)
  # nodes = del.tri$trlist[,1:3]
  # vol = (1/6)*sum(apply(nodes,1,function(loc){
  #   pp = gauge.pars[loc]
  #   ll = par.locs[loc,]
  #   coplanar.mat = do.call(rbind,lapply(c(1:3)[-1],function(i){
  #     (pp[1]*ll[1,])-(pp[i]*ll[i,])
  #   }))
  #   norm.vec = c(1,-1,1)*sapply(c(1:3),function(i){det(coplanar.mat[,-i])})
  #   return(sum(norm.vec*pp[1]*ll[1,]))
  #   # return(abs(det(cbind(rbind(0,pp*ll),1))))
  # }))
  
  del.tri = geometry::delaunayn(p=par.locs[,-3], output.options=TRUE)
  nodes = del.tri$tri
  vol2 = (1/factorial(3))*sum(apply(nodes,1,function(loc){
    pp = gauge.pars[loc]
    ll = par.locs[loc,]
    return(abs(det(ll*pp)))
  }))
  # 
  # del.tri = geometry::delaunayn(p=c(gauge.pars,0)*rbind(par.locs,0), output.options=TRUE)
  # vol3 = sum(del.tri$areas)
  # 
  # print(vol)
  # print(vol2)
  # print(vol3)
  # stop()
  
  return(vol2)
}

#################################################################################

which.adj.angles.4d = function(angles,locs){
  
  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  
  # del.tri = geometry::delaunayn(p=rbind(locs,0)[,-4], output.options=TRUE)
  del.tri = geometry::delaunayn(p=locs[,-4], output.options=TRUE)
  # w.adj.angles = tsearchn(x=rbind(locs,0)[,-4],t=del.tri$tri,xi=matrix(angles[,-4],ncol=3))
  w.adj.angles = tsearchn(x=locs[,-4],t=del.tri$tri,xi=matrix(angles[,-4],ncol=3))
  locs.idx =w.adj.angles$idx
  w.adj.angles = lapply(c(1:nrow(angles)),function(i){
    return(list(w=angles[i,],                                           # the angle
                loc.idx=del.tri$tri[locs.idx[i],],
                vertices=locs[del.tri$tri[locs.idx[i],],]))  # the enclosing vertices
  })
  
  return(w.adj.angles)
}

G.vol.4d = function(gauge.pars,par.locs){
  
  # # Method 1: determinants
  # del.tri = geometry::delaunayn(p=par.locs[,-4], output.options=TRUE)
  # nodes = del.tri$tri
  # vol = (1/factorial(4))*sum(apply(nodes,1,function(loc){
  #   pp = gauge.pars[loc]
  #   ll = par.locs[loc,]
  #   coplanar.mat = do.call(rbind,lapply(c(1:4)[-1],function(i){
  #     (pp[1]*ll[1,])-(pp[i]*ll[i,])
  #   }))
  #   norm.vec = c(1,-1,1,-1)*sapply(c(1:4),function(i){det(coplanar.mat[,-i])})
  #   return(sum(norm.vec*pp[1]*ll[1,]))
  # }))
  
  # Method 2: coordinate geometry
  del.tri = geometry::delaunayn(p=par.locs[,-4], output.options=TRUE)
  nodes = del.tri$tri
  vol2 = (1/factorial(4))*sum(apply(nodes,1,function(loc){
    pp = gauge.pars[loc]
    ll = par.locs[loc,]
    return(abs(det(ll*pp)))
  }))

  # # Method 3: using Delaunay triangulation
  # del.tri = geometry::delaunayn(p=c(gauge.pars,0)*rbind(par.locs,0), output.options=TRUE)
  # vol3 = sum(del.tri$areas)
  
  # print(vol)
  # print(vol2)
  # print(vol3)
  
  return(vol2)
}

# # A goodness-of-fit test for pwlin gauge:
# gauge.mse = function(w.mesh,pwlin.gauge,true.gauge){
#   
# }

###################################################################

# General d-dimensions

which.adj.angles = function(angles,locs){
  
  # angles are d-dimensional
  
  if("data.frame" %in% class(locs)){
    locs=as.matrix(locs)
  }
  if("data.frame" %in% class(angles)){
    angles=as.matrix(angles)
  }
  
  num.cols = dim(locs)[2]
  
  # del.tri = geometry::delaunayn(p=rbind(locs,0)[,-num.cols], output.options=TRUE)
  del.tri = geometry::delaunayn(p=locs[,-num.cols], output.options=TRUE)
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

# Given d-collumn exceedance angles, get reference angles such that each has 
# sufficient data to estmate
get.par.logs = function(data){

  
  
  return(res)
}

###################################################################

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

cart2pol.L1Rd = function(x){
  if(dim(x)[2]!=2){
    stop("input needs to be a 2 column matrix")
  }
  r<-apply(x,1,function(xx) sum(abs(xx)))
  w.cart<-x/r
  eps.val = ifelse(w.cart[,2]>=0,1,-1)
  w.pol=(1-w.cart[,1])*sign(w.cart[,2])
  # w.pol[r==0] = 0
  return(list(r=r,w=w.pol))
}

cos1 = function(q){
  q = q %% 4
  q[q>2] = q[q>2]-4
  return(1-abs(q))
}
sin1 = function(q){
  return(cos1(q-1))
}
pol2cart.L1Rd = function(r=NULL,w){
  w1 = cos1(w)
  w2 = cos1(w-1)
  if(is.null(r)){
    return(cbind(w1,w2))
  } else {
    return(r*cbind(w1,w2))
  }
}
