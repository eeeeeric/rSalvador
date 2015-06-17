# rSalvador: An R tool for the Luria-Delbruck fluctuation assay
# Qi Zheng, Department of epidemiology and Biostatistics
# Texas A&M School of Public Health
# Version 1.0: April 20, 2014
# Version 1.1: April 19, 2015
# Version 1.2: June 2, 2015

# For non-essential functions, June 5, 2015

# ----------- plating likelihood under Mandelbrot-Koch, June 5, 2015  ----------

plot.likelihood.MK=function(data,w=1,m.low=-1,m.up=-1, init.m=0, plot.pts=30,lik.col='black',
mle.col='red', title='', x.lab='Value of m', y.lab='Log-likelihood',show.secant=TRUE) {

if (title=='') {tit=paste('Log-likelihood function for', deparse(substitute(data)),
  '( w=',toString(w),')')}
else tit=title

mle=newton.MK(data,w,init.m=init.m)

if (m.low==-1) {m0=mle*0.7}
 else {m0=m.low}

if (m.up==-1) {m1=mle*1.3}
 else {m1=m.up}

max.m=log.likelihood.MK(data,mle,w)

m.pts=seq(m0,m1,length.out=plot.pts)

likely=sapply(m.pts,log.likelihood.MK,data=data,w=w)

plot(m.pts,likely,type='l',col=lik.col,main=tit,xlab=x.lab,ylab=y.lab)

abline(v=mle, col=mle.col)

if (show.secant) abline(h=max.m, lty=2)

grid()

}

