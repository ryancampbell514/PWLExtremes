require(evd)

# normalized densities
h.log = function(w,dep){
  (1/2)*evd::hbvevd(x = w, dep = dep, model="log")
}
h.neglog = function(w,dep){
  (1/2)*evd::hbvevd(x = w, dep = dep, model="neglog")
}
h.bilog = function(w,dep){
  (1/2)*evd::hbvevd(x = w, alpha=dep[1], beta=dep[2], model="bilog")
}
h.negbilog = function(w,dep){
  (1/2)*evd::hbvevd(x = w, alpha = dep[1], beta=dep[2], model="negbilog")
}
h.ct = function(w,dep){
  (1/2)*evd::hbvevd(x = w, alpha = dep[1], beta=dep[2], model="ct")
}

# spectral (probability) measures
H.log.2d = function(w,dep){
  sapply(w, function(ww) integrate(h.log,lower=0,upper=ww,dep=dep)$value)
}
H.neglog.2d = function(w,dep){
  sapply(w, function(ww) integrate(h.neglog,lower=0,upper=ww,dep=dep)$value)
}
H.bilog.2d = function(w,dep){
  sapply(w, function(ww) integrate(h.bilog,lower=0,upper=ww,dep=dep)$value)
}
H.negbilog.2d = function(w,dep){
  sapply(w, function(ww) integrate(h.negbilog,lower=0,upper=ww,dep=dep)$value)
}
H.ct.2d = function(w,dep){
  sapply(w, function(ww) integrate(h.ct,lower=0,upper=ww,dep=dep)$value)
} 

# negative log-likelihoods
nll.log = function(par,wexc){
  # dep - length 1 between 0 and 1
  if(par > 1 | par <= 0){
    return(Inf)
  } else {
    return(-sum(log(h.log(w = wexc, dep = par))))
  }
}
nll.neglog = function(par,wexc){
  # dep - length 1 between 0 and 1
  if(par<=0){
    return(Inf)
  } else {
    -sum(log(h.neglog(w = wexc, dep = par)))
  }
}
nll.bilog = function(par,wexc){
  # dep - length 2 between 0 and 1
  if(any(par<0)|any(par>1)){
    return(Inf)
  } else {
    -sum(log(h.bilog(w = wexc, dep = par)))
  }
}
nll.negbilog = function(par,wexc){
  # dep - length 2 between 0 and 1
  if(any(par<0)){
    return(Inf)
  } else {
    -sum(log(h.negbilog(w = wexc, dep = par)))
  }
}
nll.ct = function(par,wexc){
  # dep - length 2 between 0 and 1
  if(any(par<0)){
    return(Inf)
  } else {
    -sum(log(h.ct(w = wexc, dep = par)))
  }
} 
