# require(geometry)
# 
# pw_lin_gauge_2d_exp = function(w,par){
#   # w   -> a scalar on the unit simplex
#   # par -> radii of length=length(mesh)
#   
#   # notice here that the mesh is uniquely defined by the number of parameters
#   # Could make the choice to define reference angles
#   
#   if(length(w)==2){
#     w = w[1]  # want simplex, not Cartesian
#   }
#   
#   mesh.len = length(par)
#   w.mesh = seq(0,1,length.out=mesh.len)
#   
#   if(w %in% w.mesh){
#     ind = which(w == w.mesh)
#   } else {
#     ind = sapply(1:mesh.len, function(ii){
#       w >= w.mesh[ii] & w < w.mesh[ii+1]
#     })
#     ind = cbind(1:mesh.len,ind)
#     ind = ind[which(ind[,2]==1),1]
#   }
#   
#   if(ind==length(par)){
#     r1 = par[ind-1]
#     r2 = par[ind]
#     w.low = w.mesh[ind-1]
#     w.up = w.mesh[ind]
#   } else {
#     r1 = par[ind]
#     r2 = par[ind+1]
#     w.low = w.mesh[ind]
#     w.up = w.mesh[ind+1]
#   }
#   
#   dx = r2*w.up - r1*w.low
#   dy = r2*(1-w.up) - r1*(1-w.low)
#   a = dy/dx
#   b = r1*(1-w.low - a*w.low)
#   if(dx==0){
#     return(w/(r1*w.low))
#   } 
#   else{
#     return((1-w-(a*w))/b)
#   }
# }
# 
# pw_lin_gauge_2d_lap = function(w,par){
# 
#   # want to be in [0,2*pi)
#   if(length(w)==2){
#     w = cart2pol(matrix(w,nrow=1))[,1]
#   }
#   if(w<0){
#     w = w + 2*pi
#   } else if(w>=2*pi){
#     w = w - 2*pi
#   }
#   
#   mesh.len = length(par)
#   w.mesh = seq(0,2*pi*(1-(1/mesh.len)),length.out=mesh.len)
#   w.mesh = c(w.mesh,2*pi)
#   
#   if(w %in% w.mesh){
#     ind = which(w == w.mesh)
#   } else {
#     ind = sapply(1:mesh.len, function(ii){
#       w >= w.mesh[ii] & w < w.mesh[ii+1]
#     })
#     ind = cbind(1:mesh.len,ind)
#     ind = ind[which(ind[,2]==1),1]
#   }
#   
#   if(ind==mesh.len){
#     r1 = par[ind]
#     r2 = par[1]
#     w.low = w.mesh[ind]
#     w.up = w.mesh[ind+1]
#   } else {
#     r1 = par[ind]
#     r2 = par[ind+1]
#     w.low = w.mesh[ind]
#     w.up = w.mesh[ind+1]
#   }
#   
#   dx = r2*cos(w.up) - r1*cos(w.low)
#   dy = r2*sin(w.up) - r1*sin(w.low)
#   a = dy/dx
#   b = r1*(sin(w.low) - a*cos(w.low))
#   
#   if(dx==0){
#     return(cos(w)/(r1*cos(w.low)))
#   } 
#   else{
#     return((sin(w)-a*cos(w))/b)
#   }
# }
# 
# # www=seq(0,2*pi,by=0.01)
# # par=model.fit.pwl$mle[2:length(model.fit.pwl$mle)]#runif(8)
# # plot(cbind(cos(www),sin(www))/sapply(www,pw_lin_gauge_2d_lap,par),type="l")
# 
# which.region = function(w,nnodes,norm=c("1","2")){
#   if(norm=="1"){
#     if(length(w)==2){
#       w = w[1]  # want simplex, not Cartesian
#     }
#     mesh.len = length(par)
#     w.mesh = seq(0,1,length.out=mesh.len)
#     if(w %in% w.mesh){
#       ind = which(w == w.mesh)
#     } else {
#       ind = sapply(1:mesh.len, function(ii){
#         w >= w.mesh[ii] & w < w.mesh[ii+1]
#       })
#       ind = cbind(1:mesh.len,ind)
#       ind = ind[which(ind[,2]==1),1]
#     }
#     return(ind)
#   } else if(norm=="2"){
#     if(length(w)==2){
#       w = cart2pol(matrix(w,nrow=1))[,1]
#       if(w<0){
#         w = w + 2*pi
#       }
#     }
#     mesh.len = length(par)
#     w.mesh = seq(0,2*pi*(1-(1/mesh.len)),length.out=mesh.len)
#     w.mesh = c(w.mesh,2*pi)
#     if(w %in% w.mesh){
#       ind = which(w == w.mesh)
#     } else {
#       ind = sapply(1:mesh.len, function(ii){
#         w >= w.mesh[ii] & w < w.mesh[ii+1]
#       })
#       ind = cbind(1:mesh.len,ind)
#       ind = ind[which(ind[,2]==1),1]
#     }
#     return(ind)
#   }
# }
# 
# pwlin_gauge2d_fitting = function(w,par,ind,norm=c("1","2")){
#   if(norm=="1"){
#     if(ind==length(par)){
#       r1 = par[ind-1]
#       r2 = par[ind]
#       w.low = w.mesh[ind-1]
#       w.up = w.mesh[ind]
#     } else {
#       r1 = par[ind]
#       r2 = par[ind+1]
#       w.low = w.mesh[ind]
#       w.up = w.mesh[ind+1]
#     }
#     
#     dx = r2*w.up - r1*w.low
#     dy = r2*(1-w.up) - r1*(1-w.low)
#     a = dy/dx
#     b = r1*(1-w.low - a*w.low)
#     if(dx==0){
#       return(w/(r1*w.low))
#     } 
#     else{
#       return((1-w-(a*w))/b)
#     }
#   } else if(norm=="2"){
#     if(ind==mesh.len){
#       r1 = par[ind]
#       r2 = par[1]
#       w.low = w.mesh[ind]
#       w.up = w.mesh[ind+1]
#     } else {
#       r1 = par[ind]
#       r2 = par[ind+1]
#       w.low = w.mesh[ind]
#       w.up = w.mesh[ind+1]
#     }
#     
#     dx = r2*cos(w.up) - r1*cos(w.low)
#     dy = r2*sin(w.up) - r1*sin(w.low)
#     a = dy/dx
#     b = r1*(sin(w.low) - a*cos(w.low))
#     
#     if(dx==0){
#       return(cos(w)/(r1*cos(w.low)))
#     } 
#     else{
#       return((sin(w)-a*cos(w))/b)
#     }
#   }
# }
# 
# pwlin_gauge2d = function(w,par,norm=c("1","2")){
#   if(norm=="1"){
#     return(pw_lin_gauge_2d_exp(w,par))
#   } else if(norm=="2"){
#     return(pw_lin_gauge_2d_lap(w,par))
#   }
# }

#############################################################################
################################ d=2 ########################################
#############################################################################

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

gfun.simple.d2pwlin.L2 = function(w,par,par.locs){
  # xyz      = 
  # par      = vector of length 2
  # par.locs = vector of length 2
  
  r1 = par[1]
  r2 = par[2]
  
  w.sort = sort(par.locs)
  w.low = w.sort[1]
  w.up = w.sort[2]
  
  dx = r2*cos(w.up) - r1*cos(w.low)
  dy = r2*sin(w.up) - r1*sin(w.low)
  a = dy/dx
  b = r1*(sin(w.low) - a*cos(w.low))
  if(dx==0){
    return(cos(w)/(r1*cos(w.low)))
  } 
  else{
    return((sin(w)-(a*cos(w)))/b)
  }
}

gfun.simple.d2pwlin.L1.Lap = function(w,par,par.locs){
  # xyz      = 
  # par      = vector of length 2
  # par.locs = vector of length 2
  
  r1 = par[1]
  r2 = par[2]
  
  # w.sort = par.locs#sort(par.locs)
  # w.low = w.sort[1]
  # w.up = w.sort[2]
  w.low = par.locs[1]
  w.up = par.locs[2]
  
  dx = r2*cos1(w.up) - r1*cos1(w.low)
  dy = r2*sin1(w.up) - r1*sin1(w.low)
  a = dy/dx
  b = r1*(sin1(w.low) - a*cos1(w.low))
  if(dx==0){
    return(cos1(w)/(r1*cos1(w.low)))
  } 
  else{
    return((sin1(w)-(a*cos1(w)))/b)
  }
}

pwlin.g.vals.2d = function(w.adj.angles,par,par.locs,norm=NULL,marg="pos"){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    which.angles = lst$idx.locs
    if(marg=="pos"){
      return(gfun.simple.d2pwlin.L1(w=lst$w,par=par[which.angles],par.locs=par.locs[which.angles]))
    } else if(marg=="Rd"){
      return(gfun.simple.d2pwlin.L1.Lap(w=lst$w,par=par[which.angles],par.locs=par.locs[which.angles]))
    }
  })
}

gfun.2d = function(x, par, ref.angles, norm=NULL, marg="pos"){
  # x-> matrix of 2 columns
  if(is.null(dim(x)) & length(x)==1 & marg=="pos"){
    x = cbind(x,1-x)
  } else if(is.null(dim(x)) & marg=="Rd"){
    x = pol2cart.L1Rd(w=x)
  } 
  if(is.null(dim(x)) & length(x)>1){
    x = matrix(x,nrow=1)
  }
  
  # if(norm=="1"){
  #   rad.inp = apply(x,1,sum)
  #   angle.inp = x[,1] / rad.inp
  # } else if(norm=="2"){
  #   pol.tranf = cart2pol(x)
  #   rad.inp = pol.tranf[,2]
  #   angle.inp = pol.tranf[,1]
  #   angle.inp = ifelse(angle.inp<0, angle.inp+(2*pi), angle.inp)
  # }
  if(marg=="pos"){
    rad.inp = apply(x,1,sum)
    angle.inp = x[,1] / rad.inp
  } else if(marg=="Rd"){
    pol.tranf = cart2pol.L1Rd(x)
    rad.inp = pol.tranf$r
    angle.inp = pol.tranf$w
  }
  
  w.adj.angles = which.adj.angles.2d(angles=angle.inp, locs=ref.angles, norm=NULL, marg=marg)
  
  gval = rad.inp*pwlin.g.vals.2d(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles,norm=NULL,marg=marg)
  
  return(gval)
}

# pwlin.g.vals.2d.lap = function(w.adj.angles,par,par.locs){
#   # w.adj.angles  -> output of which.adj.angles()
#   # par           -> parameters at reference angles
#   # par.locs      -> 3-column matrix of reference angles.
#   
#   sapply(w.adj.angles, function(lst){
#     which.angles = lst$idx.locs
#     return(gfun.simple.d2pwlin.lap(w=lst$w,par=par[which.angles],par.locs=par.locs[which.angles]))
#   })
# }



#############################################################################
################################ d=3 ########################################
#############################################################################

# given 3 parameters (par) located at 3 reference angles (par.locs),
# return the gauge function value at the point xyz
gfun.simple.d3pwlin = function(xyz,par,par.locs=diag(3)){
  # par      = in R^3_+
  # par.locs = angles (in the simplex) where the parameters are location
  norm.vec = pracma::cross(par[1]*par.locs[1,] - par[2]*par.locs[2,],
                           par[1]*par.locs[1,] - par[3]*par.locs[3,])  # TODO: replace wit determinant
  sum(norm.vec * xyz) / sum(norm.vec * par[1] * par.locs[1,])
}


# evaluating the gauge at a set of angles, returns a vector
pwlin.g.vals.3d = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    angle = lst$w
    which.angles = lst$idx.locs
    if(any(is.na(angle)) | any(angle<0) | all(which.angles==0) | any(is.na(which.angles))){  # if angle is not in the positive orthant OR if angle lies on line
      return(NA)
    } else {
      return(gfun.simple.d3pwlin(xyz=angle,par=par[which.angles],par.locs=par.locs[which.angles,]))
    }
  })
}

gfun.3d = function(x, par, ref.angles){
  # x-> matrix of 3 columns, a bunch of angles
  # par -> a vector of distances
  # ref.angles -> a matrix of 3 columns, the reference angles
  
  if(is.null(dim(x))){
    x = matrix(x,nrow=1)
  }
  
  # del.tri = geometry::delaunayn(p=rbind(ref.angles,0)[,-4], output.options=TRUE)
  rad.inp = apply(x,1,sum)
  angle.inp = x / rad.inp
  
  # if("data.frame" %in% class(angle.inp)){
  #   angle.inp = as.matrix(angle.inp)
  # }
  # if("data.frame" %in% class(ref.angles)){
  #   ref.angles = as.matrix(ref.angles)
  # }
  
  w.adj.angles = which.adj.angles.3d(angles=angle.inp, locs=ref.angles)
  
  gval = rad.inp*pwlin.g.vals.3d(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)
  
  return(gval)
}

pwlin.g.vals.3d.lap = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    which.angles = lst$idx.locs
    return(gfun.simple.d3pwlin(xyz=lst$w,par=par[which.angles],par.locs=par.locs[which.angles,]))
  })
}

############################################################################

gfun.simple.d4pwlin = function(xyz,par,par.locs=diag(4)){
  # par      = in R^3_+
  # par.locs = angles (in the simplex) where the parameters are location
  coplanar.mat = do.call(rbind,lapply(c(1:4)[-1],function(i){
    (par[1]*par.locs[1,])-(par[i]*par.locs[i,])
  }))
  norm.vec = c(1,-1,1,-1)*sapply(c(1:4),function(i){det(coplanar.mat[,-i])})
  # TODO: replace with determinant
  # TODO: compare with pracma::crossn
  sum(norm.vec * xyz) / sum(norm.vec * par[1] * par.locs[1,])
}

# evaluating the gauge at a set of angles, returns a vector
pwlin.g.vals.4d = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    which.angles = lst$loc.idx
    if(any(is.na(lst$w)) | all(which.angles==0) | any(lst$w<0) | any(is.na(which.angles))){  # if angle is not in the positive orthant OR if angle lies on line
      return(NA)
    } else {
      return(gfun.simple.d4pwlin(xyz=lst$w,par=par[which.angles],par.locs=par.locs[which.angles,]))
    }
  })
}

# now evaluate on Cartesian points
gfun.4d = function(x, par, ref.angles){
  # x-> matrix of 4 columns
  if(is.null(dim(x))){
    x = matrix(x,nrow=1)
  }
  
  # del.tri = geometry::delaunayn(p=rbind(ref.angles,0)[,-4], output.options=TRUE)
  rad.inp = apply(x,1,sum)
  angle.inp = x / rad.inp
  is.valid = apply(angle.inp,1,function(vec) all(!is.na(vec)))
  w.adj.angles = replicate(nrow(angle.inp), 
                           list(w=rep(NA,4),  
                                loc.idx=rep(NA,4),
                                vertices=matrix(NA,4,4)), FALSE)
  inn = as.matrix(angle.inp[is.valid,],ncol=4,byrow=T)
  if(ncol(inn)==1){
    inn = t(inn)
  }
  w.adj.angles[is.valid] = which.adj.angles.4d(angles=inn, locs=ref.angles)
  
  gval = rad.inp*pwlin.g.vals.4d(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)
  return(gval)
}

############################################################################

# General d-dimensions

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
  
  # # TODO: compare with pracma::crossn
  # print(norm.vec)
  # 
  # print(par)
  # print(par.locs)
  # print(par*par.locs)
  # print(pracma::crossn((par*par.locs)[-1,]))
  # stop()
  
  
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

# now evaluate on Cartesian points
gfun.pwl = function(x, par, ref.angles){
  # x-> matrix of d columns
  if(is.null(dim(x))){
    x = matrix(x,nrow=1)
  }
  
  num.cols = dim(x)[2]
  
  # del.tri = geometry::delaunayn(p=rbind(ref.angles,0)[,-4], output.options=TRUE)
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

############################################################################
## Parametric Gauges in Laplace margins for general d-dimensions (not in geometricMVE)
gauge_gaussian_lap = function(vec,Sigma){
  Q = solve(Sigma)
  dd = length(vec)
  sgn = sign(vec)
  vec.mat = matrix(vec,ncol=1)
  sgn.mat = matrix(sgn,ncol=1)
  t(sgn.mat * sqrt(abs(vec.mat))) %*% Q %*% (sgn.mat * sqrt(abs(vec.mat)))
}
gauge_logistic_lap = function(vec,dep.par){
  dd = length(vec)
  is.pos = vec>0
  is.neg = vec<0
  is.zero = vec==0
  if(sum(is.pos)==dd){
    return((1/dep.par)*sum(vec) + (1-(dd/dep.par))*min(vec))
  } else if(sum(is.neg)==dd){
    return((sum((-vec)^(1/dep.par)))^dep.par)
  } else {
    return((1/dep.par)*sum(vec[is.pos]) + (sum((-vec[is.neg])^(1/dep.par)))^dep.par)
  }
}
##########################################################################
gauss.alog.mixgauge = function(xyz,par){
  # gauge for d=3 simulation study
  # par -> length 5
  if(length(par) != 5){
    stop("parameters not of correct length")
  }
  min(geometricMVE::gauge_gaussian3d(xyz,par=par[1:3]),
      geometricMVE::gauge_rvad_all3d(xyz,par=par[4:5],
                                     theta1=0,theta2=0,theta3=0,theta13=0,theta23=0))
}

gauss.mixgauge = function(xyz,par){
  # gauge for d=3 simulation study
  # par -> length 5
  if(length(par) != 6){
    stop("parameters not of correct length")
  }
  min(geometricMVE::gauge_gaussian3d(xyz,par=par[1:3]),
      geometricMVE::gauge_gaussian3d(xyz,par=par[4:6]))
}
