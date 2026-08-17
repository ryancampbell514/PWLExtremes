# d=3,4,5 Delauney triangulation on air pollution
rm(list = ls())

library(geometry)
library(pracma)
library(rgl)
library(openair)
library(lubridate)
library(dplyr)

setwd("~/Dropbox/phd_research/pw_lin_gauge/pollution_example/")

# fn.dir = "~/Dropbox/phd_research/pw_lin_gauge/Functions"
# invisible(sapply(file.path(fn.dir,list.files(fn.dir)),source))

# the perspective  at which we want to view rgl figures
usermat = matrix(c(-0.91820633,0.3960208,-0.008041704,0,
                   -0.09606559,-0.2029479,0.974465668,0,
                   0.38427672,0.8955331,0.224392071,0,
                   0,0,0,1)
                 ,4,4,byrow=T) 

# View(importMeta(source = "aurn", all = TRUE))

##############################################################################

# import the data for site "kc1", North Kensington, London
ds = importAURN(site = "kc1", year = 1960:2024) # "kc1" -> 1,3,7,8
head(ds)
range(ds$date)
nrow(ds)

ds.reduced = ds %>%
  group_by(Day = as.Date(date)) %>%
  dplyr::summarise(co = max(co,na.rm=T), no2=max(no2,na.rm=T), pm10=max(pm10,na.rm=T), no=max(no,na.rm=T))
ds.reduced = ds.reduced[!apply(ds.reduced,1,function(vec) any(is.na(vec))),]
ds.reduced = ds.reduced[apply(ds.reduced,1,function(vec) all(is.finite(as.numeric(vec[-1])))),]
ds.df = as.data.frame(ds.reduced[,-1])
par(mfrow=c(2,2))
for(col in 1:4){acf(ds.df[,col])}

head(ds.df)

# transform to Exp(1) margins
thresh.vals = sapply(ds.df,function(x) unname(quantile(x,0.95)))
gpd.pars = data.frame(mapply(function(x,threshold){
  out = evd::fpot(x=x,threshold=threshold)$estimate
  return(c(loc=threshold,out))
}, x=ds.df,threshold=thresh.vals))
emp.cdf<-function(x,gpd.pars){
  is.exc = x>gpd.pars[1]
  x.unif = rank(x)/(length(x)+1)
  x.unif[is.exc] = 1-mean(is.exc)*(sapply(1+gpd.pars[3]*((x[is.exc]-gpd.pars[1])/gpd.pars[2]) ,function(xx) max(xx,0)))^(-1/gpd.pars[3])
  return(x.unif)
}
ds.exp = -log(1-mapply(function(x,pars) emp.cdf(x,pars),
                        x=ds.df,pars=gpd.pars))
head(ds.exp)


# collect all the faces in a list
pollutants = colnames(ds.exp)
all.faces = lapply(1:ncol(ds.exp), function(x) combn(ncol(ds.exp),x))
all.faces[[1]] = matrix(all.faces[[1]],nrow=1)
all.faces = lapply(all.faces,function(x) split(x, rep(1:ncol(x), each = nrow(x))))
all.faces = unlist(all.faces, recursive=F, use.names=F)

# compute empirical chi(u) for u=0.95
chi.u.vals = sapply(all.faces,function(idx,u=0.98){
  if(length(idx)==1){
    return(NA)
  } else {
    dat = pexp(ds.exp[,idx])
    (1/(1-u))*mean(apply(dat,1,function(vec) all(vec>u)))
  }
})
all.faces[chi.u.vals>0]           # AD
all.faces[which.max(chi.u.vals)]  # max AD
faces[chi.u.vals==0]          # AI
## It seems like we only have strong AD (large chi values) or negative dependence with O3

plot3d(ds.exp[,-4])  # add 8 later on

# create and save the d=3,4 datasets
dim(ds.exp)
ds.exp.3d = ds.exp[,-4]
ds.exp.4d = ds.exp
save(ds.exp.3d,file="ds_3d_urban.Rdata")
save(ds.exp.4d,file="ds_4d_urban.Rdata")
par(pty="s")
plot3d(ds.exp.3d)

##############################################################################
##############################################################################

# find something more rural, Ladybower, East Midlands
ds = importAURN(site = "LB", year = 1960:2024) # "kc1" -> 1,3,7,8
head(ds)
range(ds$date)
nrow(ds)
is.winter = month(ds$date) >= 10 |  month(ds$date) <= 4
ds.winter = ds[is.winter,]
ds = ds[apply(ds,1,function(vec) any(!is.na(vec[-c(1:4)]))),]

par(mfrow=c(3,3))
for(col in 5:ncol(ds)){acf(ds[,col])}
ds = ds[seq(1,nrow(ds),by=24),]  # 24 hour - daily
par(mfrow=c(3,3))
for(col in 5:ncol(ds)){acf(ds[,col])}
dim(ds)
head(ds)
ds = as.data.frame(ds[,c(5,6,7,8,9,12)])
head(ds)

thresh.vals = sapply(ds,function(x) unname(quantile(x,0.95)))
gpd.pars = data.frame(mapply(function(x,threshold){
  out = evd::fpot(x=x,threshold=threshold)$estimate
  return(c(loc=threshold,out))
}, x=ds,threshold=thresh.vals))
emp.cdf<-function(x,gpd.pars){
  is.exc = x>gpd.pars[1]
  x.unif = rank(x)/(length(x)+1)
  x.unif[is.exc] = 1-mean(is.exc)*(sapply(1+gpd.pars[3]*((x[is.exc]-gpd.pars[1])/gpd.pars[2]) ,function(xx) max(xx,0)))^(-1/gpd.pars[3])
  return(x.unif)
}
ds.exp = -log(1-mapply(function(x,pars) emp.cdf(x,pars),
                       x=ds,pars=gpd.pars))

pollutants = colnames(ds.exp)
all.faces = lapply(1:ncol(ds.exp), function(x) combn(ncol(ds.exp),x))
all.faces[[1]] = matrix(all.faces[[1]],nrow=1)
all.faces = lapply(all.faces,function(x) split(x, rep(1:ncol(x), each = nrow(x))))
all.faces = unlist(all.faces, recursive=F, use.names=F)
faces = all.faces[1:56]

chi.u.vals = sapply(faces,function(idx,u=0.98){
  if(length(idx)==1){
    return(NA)
  } else {
    dat = pexp(ds.exp[,idx])
    (1/(1-u))*mean(apply(dat,1,function(vec) all(vec>u)))
  }
})
faces[chi.u.vals>0]           # AD
faces[which.max(chi.u.vals)]  # max AD
faces[which.min(chi.u.vals)]  # min AD
faces[chi.u.vals==0]          # AI

# 2,3,5,6
plot3d(ds.exp[,c(2,3,5)])

# create and save the d=3,4 datasets
head(ds.exp)
dim(ds.exp)
ds.exp.3d = ds.exp[,c(2,3,5)]
ds.exp.4d = ds.exp[,c(2,3,5,6)]
save(ds.exp.3d,file="ds_3d_rural.Rdata")
save(ds.exp.4d,file="ds_4d_rural.Rdata")
par(pty="s")
plot3d(ds.exp.3d)

