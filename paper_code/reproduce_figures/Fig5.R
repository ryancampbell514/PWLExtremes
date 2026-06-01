rm(list=ls())
library(geometricMVE)

load("~/Dropbox/phd_research/pw_lin_gauge/angle_density/d3_mix_fW_ex_pen_lesslocs.Rdata")

w.mesh = qr$wpts
mle = fW.fit.pen$fW.mle
g.vals.fW = gfun.pwl(x=w.mesh,par=mle,#W.pars.34[best.mod.idx,],
                     ref.angles=par.locs)
G.vol.W = G.vol.3d(mle,#W.pars.34[best.mod.idx,],
                   par.locs)
d = data.frame(x=w.mesh[,1], y=w.mesh[,2], z=(1/(3*G.vol.W))*(g.vals.fW^(-3)))

pdf("figures/d3_mix_fWfit.pdf",width=5,height=5)
par(pty="s",mar = c(4, 4.5, 0.1, 0.1),cex.lab=1.3)
lattice::levelplot(z~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))
trellis.focus("panel", 1, 1, highlight=FALSE)
lpoints(par.locs[,1], par.locs[,2], pch=20, col="maroon1", cex=3)
trellis.unfocus()
dev.off()

max.clip = max((1/(3*G.vol.W))*(g.vals.fW^(-3)),na.rm=T)
wexc.samp.density = squash::hist2(x=wexc[,1],y=wexc[,2],plot=F)
xy = expand.grid(wexc.samp.density$x[-length(wexc.samp.density$x)],
                 wexc.samp.density$y[-length(wexc.samp.density$y)])
d = data.frame(x=xy[,1], y=xy[,2], 
               z=as.numeric(wexc.samp.density$z))
d$dens = d$z / (sum( d$z, na.rm=T)*(0.02^2))
d$dens.clip = ifelse(d$dens>max.clip,max.clip,d$dens)
pdf("figures/d3_mix_Wexc_samp.pdf",width=5,height=5)
par(pty="s",mar = c(4, 4.5, 0.1, 0.1),cex.lab=1.3)
lattice::levelplot(dens.clip~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))
dev.off()


max.clip = max((1/(3*G.vol.W))*(g.vals.fW^(-3)),na.rm=T)
wexc.samp.density = squash::hist2(x=w[,1],y=w[,2],plot=F)
xy = expand.grid(wexc.samp.density$x[-length(wexc.samp.density$x)],
                 wexc.samp.density$y[-length(wexc.samp.density$y)])
d = data.frame(x=xy[,1], y=xy[,2], 
               z=as.numeric(wexc.samp.density$z))
d$dens = d$z / (sum( d$z, na.rm=T)*(0.02^2))
d$dens.clip = ifelse(d$dens>max.clip,max.clip,d$dens)
pdf("figures/d3_mix_W_samp.pdf",width=5,height=5)
par(pty="s",mar = c(4, 4.5, 0.1, 0.1),cex.lab=1.3)
lattice::levelplot(dens.clip~x*y, data=d,main=NA,xlab=expression(w[1]),ylab=expression(w[2]))
dev.off()
