# rSalvador: An R tool for the Luria-Delbruck fluctuation assay
# Qi Zheng, Department of epidemiology and Biostatistics
# Texas A&M School of Public Health
# Version 1.0: April 20, 2014
# Version 1.1: April 19, 2015
# Version 1.2: June 2, 2015


# ----------- July 9, 2013, first successful calling C from R for Salvador
# ----------- Version 0.2 started December 22, 2013
library(gdata)
library(hypergeo)

# --------------- user interface, July 17, 2013

# ------  excel read, July 17, 2013

import.excel.data=function(filename,col=1) {

x=read.xls(filename,header=T)
y=x[[col]]
n=length(y)
mutant.name=names(x)[col]
message('')
message(paste('Experiment is labeled ',trim(toString(mutant.name)),'-- it consists of ',
     toString(n), ' cultures. Observations are'))
message('')
print(y)

return(y)
}

# --------- text read, July 17, 2013
import.text.data=function(filename,jump=1,col=1) {

x=read.table(filename,skip=jump)

y=as.numeric(unlist(x[,col]))

n=length(y)

message('')
message(paste('Experiment  consists of ', toString(n), ' cultures. Observations are'))
message('')
print(y)

return(y)
}

# -------------- plain text write, July 13, 2013

export.text.data=function(filename,data) {

f=file(filename,'w')

myheader=paste('Fluctuation Assay by rSalvador:', toString(date()) )

writeLines(myheader,f)

write.table(data,f,row.names=F,col.names=F)

close(f)

message('')
message(paste('Data saved to file ',toString(filename),' ...'))
message('')
}



# ------------- simuMutants, July 10, 2013

simu.mutants=function(mu,b1,b2,N0,Nt) {

T=log(Nt/N0)/b1

mean.mutation=mu*N0*(exp(b1*T)-1)/b1

n.mutation=rpois(1,mean.mutation)

mutation.epoch=log(1+runif(n.mutation)*(exp(b1*T)-1))/b1

offspring=lapply(mutation.epoch,function(x) {1+rgeom(1,exp(-b2*(T-x)))} )

offspring=sum(unlist(offspring))

return(offspring)

}

# ------------------------ simulation of all sister cultures in an experiment

simu.cultures=function(n,mu,b1,b2,N0,Nt) {

cultures=rep(-1,n)

for (i in 1:n) {

cultures[i]=simu.mutants(mu,b1,b2,N0,Nt) }

return(cultures)

}



# --------- the p0 method ------------

LD.p0.est=function(data) {

p0=sum(data==0)/length(data)

if (p0>0) {return (-log(p0)) }

else {message('The P0 method is not applicable.')}

}



# --------------- Jones Median estimator 

jones.median.est=function(data) {

r=median(data)

est=(r-0.693147)/( log(r)+0.3665)

return (est)

}



# ---------------  computing  the LD distribution ----------------
prob.LD=function(m,phi=1,k) {   ### June 3, 2015: let phi default to one.

result.p=rep(1.1,k+1)

####### if(! is.loaded('delbruck')) dyn.load('./delbruck.so')

z=.C('luria_R_wrapper',m=as.double(m),phi=as.double(phi),k=as.integer(k),result=as.double(result.p))

return(z$result)

}

# ---------------------- convolution of two sequences -----------
seq.convolute=function(x,y) {

xlen=length(x)

ylen=length(y)

result=rep(0,xlen)

u=.C('pdfconv_R_wrapper',x=as.double(x),xlen=as.integer(xlen),y=as.double(y),

    ylen=as.integer(xlen), z=as.double(result))

return(u$z)

}

# -------------------- Score and Information for LD distribution 

LD.score.info=function(m,phi,data) {

n=max(data)

p=prob.LD(m,phi,n)

h=apply(as.array(1:n), 1, function(x) phi^(x-1) *(1/x -phi/(x+1)) )

h=append(h,-1,after=0)

p1=seq.convolute(p,h)

p2=seq.convolute(p1,h)

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c(sum(score), sum(info) ) )

}

# -------------------- Score, Information and Log-likelihood for LD distribution 
# ### July 10, 2013: combining LDLoglikeScore and LDSCore$Info into one function

LD.score.info.loglikely=function(m,phi,data) {

n=max(data)

p=prob.LD(m,phi,n)

loglike=log(p[data+1])

h=apply(as.array(1:n), 1, function(x) phi^(x-1) *(1/x -phi/(x+1)) )

h=append(h,-1,after=0)

p1=seq.convolute(p,h)

p2=seq.convolute(p1,h)

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c(sum(score), sum(info), sum(loglike) ) )

}


# -------------------- LDLoglikeScore[m,phi,data] from SALVADOR 2.3, July 10, 2013

LD.loglike.score=function(m,phi,data) {

n=max(data)

p=prob.LD(m,phi,n)

loglike=log(p[data+1])

h=apply(as.array(1:n), 1, function(x) phi^(x-1) *(1/x -phi/(x+1)) )

h=append(h,-1,after=0)

p1=seq.convolute(p,h)

score=(p1/p)[data+1]

return( c(sum(loglike), sum(score)) )

}


# -------------------- the LD likelihood function, July 12, 2013

log.likelihood.LD=function(data,m,phi=1) {

n=max(data)

p=prob.LD(m,phi,n)

loglike=sum(log(p[data+1]))

return(loglike)

} # end of log.likelihood.LD




# ===================== newtonLD ===============

newton.LD=function(data,phi=1.0, tol=0.00000001, init.m=0, max.iter=30, show.iter=FALSE) {

if (init.m>0) {m0=init.m}

else { if (median(data)==0) m0=LD.p0.est(data) else m0=jones.median.est(data) }

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

for (i in 1:max.iter) {

score.info=LD.score.info(m0,phi,data)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}

# ------------ July 10, 2013, Likelihood ratio CI for LD

confint.LD=function(data,alpha=0.05,phi=1,tol=0.000001,init.m=0,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

initLow=init.lower

initUp=init.upper

mhat=newton.LD(data,phi=phi,tol=tol,init.m=init.m)

if (show.iter) {message( paste('The ML estimate of M ... ',toString(mhat) ) ) }

score.info.like=LD.score.info.loglikely(mhat,phi=phi,data=data)

score=score.info.like[1]

info=score.info.like[2]

like=score.info.like[3]

qa=qchisq(1-alpha,1)

h=sqrt(qa/info)

la=like-0.5*qa

if (show.iter) message('Iterating for lower limit ... ')

if (initLow<=0) {m0=mhat-0.5*h} else {m0=initLow}

for (i in 1:max.iter) {

 like.score=LD.loglike.score(m0,phi,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mL=m1; break} else {m0=m1}

} # end of 1st for loop


if (show.iter) message('Iterating for upper limit ... ')

if (initUp<=0) {m0=mhat+0.5*h} else {m0=initUp}

for (i in 1:max.iter) {

 like.score=LD.loglike.score(m0,phi,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mU=m1; break} else {m0=m1}

} # end of 2nd for loop

return (c(mL,mU))

}  # end of function

# -------------------- plating efficiency, July 13, 2013

etaSmall=function(e,n) {

odds=e/(1-e)

eta=rep(0,n)

eta[1]=-1-log(e)/(1-e)

for (i in 1:(n-1) ) {

eta[i+1]=-odds*eta[i]+1/i/(i+1)

}

eta=append(eta,log(e),after=0)

return(eta)

}  # end of etaSmall

etaBig=function(e,n) {

revOR=(1-e)/e

eta=rep(0,n)

eta[n]=( (1-e)/n/(n+1) ) * Re( hypergeo(1,2,n+2,1-e) )

for (i in (n-1):1 ) {

eta[i]=revOR*( 1/i/(i+1) - eta[i+1] )

}

eta=append(eta,log(e),after=0)

return(eta)

}

# -------------- computing the eta sequence -----------
etaSeq=function(e,n) {

if (e<=0.5) { eta=etaSmall(e,n) } else {eta=etaBig(e,n)}

return(eta)

}


# -------------- P0 method for plating, July 14, 2013

p0.plating=function(data,e) {

z=sum(data==0)

if (z==0) {message('P0 method is not applicable.');return()}

 else { n=length(data);  est=(1-e)*( log(z/n)/e/log(e)); return(est) }

}

# ------------- Jones median method for plating, July 14, 2013

jones.median.plating=function(data,e) {

r=median(data)

if (r==0) { message('Jones median method is not applicable.'); return()}

est=(r/e -log(2))/(log(r/e) -log(log(2)))

return(est)

}


# ---------------  plating, July 13, 2013 --------------

prob.LD.plating=function(m,e,n) {

eta.seq=etaSeq(e,n)

result.p=rep(0,n+1)

z=.C('pmfLDPlat_R_wrapper',m=as.double(m),e=as.double(e),eta=as.double(eta.seq),

   etaLen=as.integer(n+1),  prob=as.double(result.p))

return(z$prob)

}

# ----------- score and fisher info for plating, July 14, 2013 --------

LD.score.info.plating=function(m,e,eta,data) {

k=e/(1-e)

n=max(data)

p=prob.LD.plating(m,e,n)

p1=seq.convolute(eta,p)*k

p2=seq.convolute(eta,p1)*k

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c(sum(score),sum(info)) )

}




# ===================== newton.LD.plating July 14, 2013  ===============

newton.LD.plating=function(data, e, tol=0.00000001, init.m=0, max.iter=30, show.iter=FALSE) {

if (init.m>0) {m0=init.m}

else { if (median(data)==0) m0=p0.plating(data,e) else m0=jones.median.plating(data,e) }

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

eta=etaSeq(e,max(data))

for (i in 1:max.iter) {

score.info=LD.score.info.plating(m0,e,eta,data)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}


# -------------------- LDLoglikeScorePlating[m,e,eta,data] from SALVADOR 2.3, July 15, 2013

LD.loglike.score.plating=function(m,e,eta,data) {

k=e/(1-e)

n=max(data)

p=prob.LD.plating(m,e,n)

p1=seq.convolute(eta,p)*k

loglike=log(p[data+1])

score=(p1/p)[data+1]

return( c(sum(loglike), sum(score)) )

}


# -------------- combining two functions into one for efficiency, used only once, July 15, 2013

LD.score.info.loglikely.plating=function(m,e,eta,data) {

k=e/(1-e)

n=max(data)

p=prob.LD.plating(m,e,n)

p1=seq.convolute(eta,p)*k

p2=seq.convolute(eta,p1)*k

loglike=log(p[data+1])

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c( sum(score), sum(info), sum(loglike) )  )

}


# ------------ July 15, 2013, Likelihood ratio CI for LD with plating efficiency e --------

confint.LD.plating=function(data,e,alpha=0.05,tol=0.000001,init.m=0,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

initLow=init.lower

initUp=init.upper

eta=etaSeq(e,max(data))

mhat=newton.LD.plating(data,e=e,tol=tol,init.m=init.m)

if (show.iter) {message( paste('The ML estimate of M ... ',toString(mhat) ) ) }

score.info.like=LD.score.info.loglikely.plating(mhat,e=e,eta=eta,data=data)

score=score.info.like[1]

info=score.info.like[2]

like=score.info.like[3]

qa=qchisq(1-alpha,1)

h=sqrt(qa/info)

la=like-0.5*qa

if (show.iter) message('Iterating for lower limit ... ')

if (initLow<=0) {m0=mhat-0.5*h} else {m0=initLow}

for (i in 1:max.iter) {

 like.score=LD.loglike.score.plating(m0,e,eta,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mL=m1; break} else {m0=m1}

} # end of 1st for loop


if (show.iter) message('Iterating for upper limit ... ')

if (initUp<=0) {m0=mhat+0.5*h} else {m0=initUp}

for (i in 1:max.iter) {

 like.score=LD.loglike.score.plating(m0,e,eta,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mU=m1; break} else {m0=m1}

} # end of 2nd for loop

return (c(mL,mU))

}  # end of function


# -------------------- the LD likelihood function with plating, July 15, 2013

log.likelihood.LD.plating=function(data,m,e) {

n=max(data)

p=prob.LD.plating(m,e,n)

loglike=sum(log(p[data+1]))

return(loglike)

} # end of log.likelihood.LD


# ------------- functions from Falcor, July 19, 2013

falcor.confint.LD=function(data) {

C=length(data)

m=newton.LD(data)

sigma=1.225*m^(-0.315)/sqrt(C)

up.limit=log(m)+1.96*sigma*(exp(1.96*sigma))^(-0.315)

low.limit=log(m)-1.96*sigma*(exp(1.96*sigma))^(+0.315)

return( c(exp(low.limit), exp(up.limit)) )

}

falcor.LD.plating=function(data,e) {

m=newton.LD(data)

adj=(e-1)/e/log(e)

return( adj*m )

}

# modification August 6, 2013, after email from Foster

falcor.confint.LD.plating=function(data,e) {

C=length(data)

m=newton.LD(data)

adj=(e-1)/e/log(e)

m.adj=m*adj

sigma=1.225*m.adj^(-0.315)/sqrt(C)

up.limit=log(m.adj)+1.96*sigma*(exp(1.96*sigma))^(-0.315)

low.limit=log(m.adj)-1.96*sigma*(exp(1.96*sigma))^(+0.315)

return( c(exp(low.limit), exp(up.limit)) )

}


# -------- to plot the likelihood function along with the MLE of m, August 2, 2013

plot.likelihood.LD=function(data, init.m=0, m.low=-1,m.up=-1,plot.pts=30,lik.col='black',
 mle.col='red', title='', x.lab='Value of m', y.lab='Log-likelihood',show.secant=TRUE) {

mle=newton.LD(data,init.m=init.m)  ### add init.m, June 5, 2015

if (title=='') {tit=paste('Log-likelihood function for', deparse(substitute(data)))}
else tit=title


if (m.low==-1) {m0=mle*0.7}
 else {m0=m.low}

if (m.up==-1) {m1=mle*1.3}
 else {m1=m.up}

max.m=log.likelihood.LD(data,mle)

m.pts=seq(m0,m1,length.out=plot.pts)

likely=sapply(m.pts,log.likelihood.LD,data=data)

plot(m.pts,likely,type='l',col=lik.col,main=tit,xlab=x.lab,ylab=y.lab)

abline(v=mle, col=mle.col)

if (show.secant) abline(h=max.m, lty=2)

grid()  ## June 5, 2015
}

# ----------- plating efficiency ----------

plot.likelihood.LD.plating=function(data,e,m.low=-1,m.up=-1,plot.pts=30,lik.col='black',
mle.col='red', title='', x.lab='Value of m', y.lab='Log-likelihood',show.secant=TRUE) {

if (title=='') {tit=paste('Log-likelihood function for', deparse(substitute(data)))}
else tit=title

mle=newton.LD.plating(data,e)

if (m.low==-1) {m0=mle*0.7}
 else {m0=m.low}

if (m.up==-1) {m1=mle*1.3}
 else {m1=m.up}

max.m=log.likelihood.LD.plating(data,mle,e)

m.pts=seq(m0,m1,length.out=plot.pts)

likely=sapply(m.pts,log.likelihood.LD.plating,data=data,e=e)

plot(m.pts,likely,type='l',col=lik.col,main=tit,xlab=x.lab,ylab=y.lab)

abline(v=mle, col=mle.col)

if (show.secant) abline(h=max.m, lty=2)

}



# ---------------  computing  the Haldane  distribution, December 22, 2013  ----------------
prob.Haldane=function(gen, mu, n, N0=1) {

result.p=rep(1.1,n+1)

z=.C('pmfHald_R_wrapper',gen=as.integer(gen),mu=as.double(mu),n=as.integer(n), 
 N0=as.integer(N0), result=as.double(result.p))

return(z$result)

}


# ---------------  simulating the Haldane distribution, December 23, 2013  ----------------
simu.Hald=function(gen, mu) {

mut=-1

z=.C('simuHald_R_wrapper',gen=as.integer(gen),mu=as.double(mu),mut=as.integer(mut)) 

return(z$mut)

}

simu.Haldane=function(gen,mu,culture=1) {
mut=rep(-1,culture)
for(i in 1:culture) { mut[i]=simu.Hald(gen,mu) }
return(mut)
}

# ============= new to rSalvador, Dec 24, 2013 =======================

simu.Kim=function(gen, mu) {

mut=-1

z=.C('simuKimmel_R_wrapper',nGen=as.integer(gen),mutRate=as.double(mu),mutants=as.integer(mut)) 

return(z$mutants)

}

simu.Kimmel=function(gen,mu,culture=1) {
mut=rep(-1,culture)
for(i in 1:culture) { mut[i]=simu.Kim(gen,mu) }
return(mut)
}

#  ---------------------- Dec 30, 2013, with both generations reported --------- 
simu.Kim2=function(gen, mu) {

mut0=-1
mut=-1

z=.C('simuKimm2_R_wrapper',nGen=as.integer(gen),mutRate=as.double(mu),
  mutants0=as.integer(mut0), mutants=as.integer(mut)) 

return(c(z$mutants0, z$mutants))

}

simu.Kimmel2=function(gen,mu,culture=1) {
mut=matrix(-1,ncol=2,nrow=culture)
for(i in 1:culture) { mut[i,]=simu.Kim2(gen,mu) }
return(mut)
}


# ---------- Jan 1, 2014 ------------
simu.Hald2=function(gen, mu) {

mut0=-1
mut=-1

z=.C('simuHald2_R_wrapper',gen=as.integer(gen),mu=as.double(mu),mut0=as.integer(mut0),
  mut=as.integer(mut)) 

return(c(z$mut0,z$mut))

}

simu.Haldane2=function(gen,mu,culture=1) {
mut=matrix(-1,ncol=2,nrow=culture)
for(i in 1:culture) { mut[i,]=simu.Hald2(gen,mu) }
return(mut)
}



# ----------  computing derivatives of the Haldane  distribution, December 23, 2013  --------

deriv.Haldane=function(gen,mu,n,N0) {

p=rep(1.2,n+1)
p1=rep(1.2,n+1)
p2=rep(1.2,n+1)

z=.C('derivHald_R_wrapper',gen=as.integer(gen),mu=as.double(mu),n=as.integer(n), 
 N0=as.integer(N0), prob=as.double(p), prob1=as.double(p1),prob2=as.double(p2))

return(list(z$prob,z$prob1,z$prob2))

}



# ---------- reduced form of the above, just p and p1, December 24, 2013  --------

pAndP1.Haldane=function(gen,mu,n,N0) {

p=rep(1.2,n+1)
p1=rep(1.2,n+1)

z=.C('pAndP1Hald_R_wrapper',gen=as.integer(gen),mu=as.double(mu),n=as.integer(n), 
 N0=as.integer(N0), prob=as.double(p), prob1=as.double(p1))

return(list(z$prob,z$prob1))

}

# ---------- Haldane U and J , December 23, 2013 --------------

getHaldane.UandJ=function(data,g,mu,N0) {

n=max(data);

deriv=deriv.Haldane(g,mu,n,N0)

p=deriv[[1]][data+1]
p1=deriv[[2]][data+1]
p2=deriv[[3]][data+1]

U=sum(p1/p)

J=sum( (p1/p)^2 - p2/p)

return(c(U,J))

}

# ---------- Haldane likelihood and score functions, December 24, 2013 --------------

getHaldane.likAndScore=function(data,g,mu,N0) {

n=max(data);

pp1=pAndP1.Haldane(g,mu,n,N0)

p=pp1[[1]][data+1]
p1=pp1[[2]][data+1]

lik=sum(log(p))
score=sum(p1/p)

return(c(lik,score))

}





# ------------- Rossman estimator for Haldane, Adapted from SALVADOR 2.3, Dec 23, 2013
# ------- must have Nt>m, Dec 23, 2013

rossman=function(data, g, N0) {

Nt=N0*2^g
m=mean(data)
return( 2*(1.0-(1.0-m/Nt)^(1.0/g)) )
}


# -------------- newton.Haldane, Dec 23, 2013  -------------------

newton.Haldane=function(data, g, N0=1, tol=0.00000001, init.mu=0, max.iter=30, show.iter=FALSE) {

if (max(data)>=N0*2^g) {message('Maximum of data is inconsistent with model specified by g');
                        return(-1)}  # ### Jan 1, 2014

if (init.mu>0) {mu0=init.mu}

else { mu0=0.5*rossman(data,g,N0) }

if (show.iter) message( paste('iteration 0 yielding ...', toString(mu0) ) )

for (i in 1:max.iter) {

score.info=getHaldane.UandJ(data,g,mu0,N0)

mu1=mu0+score.info[1]/score.info[2]

if ( abs(mu1-mu0)/mu0 <tol ) {return (mu1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(mu1) )) ;

 mu0=mu1

}  # end of if-else

} # end of for 

return('no convergence')

}


# ------------ December 24, 2013, Likelihood ratio CI for Haldane

confint.Haldane=function(data,g,N0=1,alpha=0.05,tol=0.000001,init.mu=0,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

initLow=init.lower

initUp=init.upper

chiAlpha=qchisq(1-alpha,1)

muHat=newton.Haldane(data,g,N0,tol=tol,init.mu=init.mu)

if (show.iter) {message( paste('The ML estimate of Mu ... ',toString(muHat) ) ) }

UandJ=getHaldane.UandJ(data,g,muHat,N0)
U=UandJ[1]
J=UandJ[2]

h=0.5*sqrt(chiAlpha/J)

lik=getHaldane.likAndScore(data,g,muHat,N0)[1]

la=lik-0.5*chiAlpha

if (show.iter) message('Iterating for lower confidence limit ... ')

if (initLow<=0) {ci0=muHat-h} else {ci0=initLow}

for (i in 1:max.iter) {

 like.score=getHaldane.likAndScore(data,g,ci0,N0)

 lik=like.score[1];  U=like.score[2]

 ci1=ci0 - (lik -la)/U

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(ci1)) ) }
 # was ci0 in SALVADOR 2.3, changed to ci1, Dec 24, 2013

 if ( abs( (ci1-ci0)/ci0 )<tol ) {muLow=ci1; break} else {ci0=ci1}

} # end of 1st for loop


if (show.iter) message('Iterating for upper confidence limit ... ')

if (initUp<=0) {ci0=muHat+h} else {ci0=initUp}

for (i in 1:max.iter) {

 like.score=getHaldane.likAndScore(data,g,ci0,N0)

 lik=like.score[1];  U=like.score[2]

 ci1=ci0 - (lik -la)/U

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(ci1)) ) }
 # was ci0 in SALVADOR 2.3, changed to ci1, Dec 24, 2013

 if ( abs( (ci1-ci0)/ci0 )<tol ) {muUp=ci1; break} else {ci0=ci1}

} # end of 2nd for loop

return (c(muLow,muUp))

}  # end of function confint.Haldane


# -------------------- the Haldane likelihood function, December 25, 2013

log.likelihood.Haldane=function(data,g,mu,N0=1) {

n=max(data)

p=prob.Haldane(g,mu,n,N0)

loglike=sum(log(p[data+1]))

return(loglike)

} # end of log.likelihood.LD


# -------- to plot the likelihood function with the MLE of mu in Haldane, Decenber 25, 2013

plot.likelihood.Haldane=function(data,g, init.mu=0, mu.low=-1,mu.up=-1,plot.pts=30,lik.col='black',
 mle.col='red', title='', x.lab='Value of mu', y.lab='Log-likelihood',show.secant=TRUE) {

mle=newton.Haldane(data,g,init.mu=init.mu)

if (title=='') {tit=paste('Log-likelihood function (Haldane model) for',
  deparse(substitute(data)))}
   else tit=title


if (mu.low==-1) {mu0=mle*0.7}
 else {mu0=mu.low}

if (mu.up==-1) {mu1=mle*1.3}
 else {mu1=mu.up}

max.mu=log.likelihood.Haldane(data,g,mle)

mu.pts=seq(mu0,mu1,length.out=plot.pts)

likely=sapply(mu.pts,log.likelihood.Haldane,data=data,g=g)

plot(mu.pts,likely,type='l',col=lik.col,main=tit,xlab=x.lab,ylab=y.lab)

abline(v=mle, col=mle.col)

grid() # June 5, 2015

if (show.secant) abline(h=max.mu, lty=2)

}


### ------------ comparison of mutation rates -----------------


# Likelihood ratio test, rewritten my Mathematica code (2008) in R, April 13, 2014
# revision, April 22, 2014
# revision, April 19, 2015

# ---------- it assumes a common Nt, and hence a common phi parameter, April 19, 2015

compare.LD=function(x1,x2,init.m0=0,init.m1=0,init.m2=0,phi=1){
m1=newton.LD(x1,phi,init.m=init.m1)
m2=newton.LD(x2,phi,init.m=init.m2)
down=likely2(x1,x2,m1,m2,phi)
y=unlist(c(x1,x2))
m0=newton.LD(y,phi,init.m=init.m0)
up=likely1(y,m0,phi)
chi=-2*(up-down)
pval=1-pchisq(chi,1)
return(c(chi,pval))
}

# ----- for compare.LD

likely2=function(x1,x2,m1,m2,phi=1) {
n1=max(x1)
n2=max(x2)
p1=prob.LD(m1,phi,n1)
p1=sum(log(p1[x1 +1]))
p2=prob.LD(m2,phi,n2)
p2=sum(log(p2[x2 +1]))
return(p1+p2)
}

# ----- for compare.LD

likely1=function(x,m,phi=1) {
n=max(x)
p=prob.LD(m,phi,n)
p=sum( log(p[x+1]) )
return(p)
}

### ------ simulation under the Bartlett model, added April 19, 2015 ----

# March 14, 2015, modified April 19, 2015

simu.Bartlett=function(b1,b2,mu,N0,T,max.events=1e10,show.growth=FALSE)
{
      wild = N0
       mutants = 0 
      events = 0
      timeNow=0
      while (TRUE) {
        if (show.growth) {message(paste('Wild cell == ', toString(wild), ' mutant cell == ',
            toString(mutants) )) }
        R = (b1 + mu)*wild + b2*mutants
        timeStep = -log( runif(1,0,1) )/R
        timeNow = timeNow+timeStep
       if (timeNow > T) {return(c(wild, mutants))} else {events=events+1}
       if (runif(1,0,1) <= b1*(wild/R) ) {wild = wild + 1}  else  {mutants = mutants+1}
       if (events>max.events) {message("simulation events exceed limit ..."); return(NA)}
}   #while 
}
  

######## May 15, 2015: fitness in Luria-Delbruck

betaSeq=function(w,n) {

r=1/w

B=rep(0,n)

B[1]=1/(1+r)

for (k in 2:n) {

  B[k]=(k-1)/(k+r) * B[k-1] }

return(B)

}


# ---------------  probabilites for fitness, May 15, 2015 --------------

prob.MK=function(m,w=1,n) {

beta.seq=betaSeq(w,n)

result.p=rep(0,n+1)

########## if(! is.loaded('kochr')) dyn.load('koch.so')

z=.C('pmfMK_R_wrapper',m=as.double(m),w=as.double(w),beta=as.double(beta.seq),

   betaLen=as.integer(n+1),  prob=as.double(result.p))

return(z$prob)

}


# ===================== newton.MK May 15, 2015  ===============

newton.MK=function(data, w=1, tol=1e-8, init.m=0, max.iter=30, show.iter=FALSE) {

if (init.m>0) {m0=init.m}

# else { if (min(data)==0) m0=LD.p0.est(data) else m0=jones.median.plating(data,e=w) }

else { if (median(data)>0) m0=jones.median.plating(data,e=w) else m0=LD.p0.est(data) }

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

B=betaSeq(w,max(data))/w

h=append(B,-1,after=0)  # the h sequence

for (i in 1:max.iter) {

score.info=MK.score.info(m0,w,h,data)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}


# ----------- score and fisher info for Mandelbrot-Koch, May 15, 2015 --------

MK.score.info=function(m,w,h,data) {

n=max(data)

p=prob.MK(m,w,n)

p1=seq.convolute(h,p)

p2=seq.convolute(h,p1)

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c(sum(score),sum(info)) )

}


# -------------------- adapted from partial plating, May 15, 2015

MK.loglike.score=function(m,w,hSeq,data) {

n=max(data)

p=prob.MK(m,w,n)

p1=seq.convolute(hSeq,p)

loglike=log(p[data+1])

score=(p1/p)[data+1]

return( c(sum(loglike), sum(score)) )

}


# -------------- Mandelbrot-Koch, adapted from partial plating, May 15, 2015

MK.score.info.loglikely=function(m,w,hSeq,data) {

n=max(data)

p=prob.MK(m,w,n)

p1=seq.convolute(hSeq,p)

p2=seq.convolute(hSeq,p1)

loglike=log(p[data+1])

score=(p1/p)[data+1]

info=( (p1/p)^2 -p2/p )[data+1]

return( c( sum(score), sum(info), sum(loglike) )  )

}


# ------------ confidence interval for MK model, fitness is w, May 15, 2015 ---------

confint.MK=function(data,w=1,alpha=0.05,tol=1e-8,init.m=0,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

initLow=init.lower

initUp=init.upper

B=betaSeq(w,max(data))/w

hSeq=append(B,-1,after=0)

mhat=newton.MK(data,w=w,tol=tol,init.m=init.m)

if (show.iter) {message( paste('The ML estimate of M ... ',toString(mhat) ) ) }

score.info.like=MK.score.info.loglikely(mhat,w=w,hSeq=hSeq,data=data)

score=score.info.like[1]

info=score.info.like[2]

like=score.info.like[3]

qa=qchisq(1-alpha,1)

h=sqrt(qa/info)

la=like-0.5*qa

if (show.iter) message('Iterating for lower limit ... ')

if (initLow<=0) {m0=mhat-0.5*h} else {m0=initLow}

for (i in 1:max.iter) {

 like.score=MK.loglike.score(m0,w,hSeq,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mL=m1; break} else {m0=m1}

} # end of 1st for loop


if (show.iter) message('Iterating for upper limit ... ')

if (initUp<=0) {m0=mhat+0.5*h} else {m0=initUp}

for (i in 1:max.iter) {

 like.score=MK.loglike.score(m0,w,hSeq,data)

 like=like.score[1];  score=like.score[2]

 m1=m0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(m1)) ) }

 if ( abs( (m1-m0)/m0 )<tol ) {mU=m1; break} else {m0=m1}

} # end of 2nd for loop

return (c(mL,mU))

}  # end of function


############ Delbruck's likely average #####################
### May 21, 2015: Delbruck's likely average to be extended 

likely.average=function(data, b=1.24, classic=FALSE, init.m=0, tol=1e-9, max.iter=25,
          use.median=TRUE, show.iter=FALSE) {

if (classic) {n=length(data)} else {n=1}

if (use.median) {med=median(data)} else {med=mean(data)} 
if (med==0) {message('Median/mean is zero ...'); return( exp(-b) ) }     
if (init.m<=0) {m0=(med-log(2))/(log(med)-log(log(2)))} else {m0=init.m}     

if (show.iter) print(m0)        

for(i in  1:max.iter) {       
       m1 = m0 - (m0*log(n*m0)+b*n*m0-med)/(log(n*m0)+1+b*n);
       if (show.iter)  print(m1)
       if (abs((m0 - m1)/m0) < tol)  {return(m1)}  else {m0 = m1}  } ### end of for
        message('no convergence')  
} ### end of coulson


##########  likelihood ratio tests ######################
# June 4, 2015 -- extending LRT for comparion to the Mandelbrot-Koch model

### finding the score and observed info for the model under null hypothesis

combo.MK.score.info=function(m,X1,X2,R, w1=1,w2=1) {

B1=betaSeq(w1,max(X1))/w1

h1=append(B1,-1,after=0)

B2=betaSeq(w2,max(X2))/w2

h2=append(B2,-1,after=0)

us1=MK.score.info(m,w1,h1,X1)

us2=MK.score.info(R*m,w2,h2,X2)

u1=us1[1]

j1=us1[2]

u2=us2[1]

j2=us2[2]

u=u1+R*u2

j=j1+R*R*j2

return(c(u,j))

}

### find Mc, the combined MLE, June 4, 2015

# =====================combined newtonMK ===============

combo.newton.MK=function(X1,X2,R, w1=1.0, w2=1.0,
           tol=1e-8, init.m=0, max.iter=30, show.iter=FALSE) {

Xc=c(X1,X2)

if (init.m>0) {m0=init.m}

   else if (init.m==0) 

       {if (min(Xc)==0) m0=LD.p0.est(Xc) else m0=jones.median.est(Xc) }  ## option A

   else 

       {if (min(X1)==0) m0=LD.p0.est(X1) else m0=jones.median.est(X1) }  ## option B

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

for (i in 1:max.iter) {

score.info=combo.MK.score.info(m0,X1,X2,R,w1,w2)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}


# -------------------- the MK likelihood function, June 4, 2015

log.likelihood.MK=function(data,m,w=1) {

n=max(data)

p=prob.MK(m,w,n)

loglike=sum(log(p[data+1]))

return(loglike)

} # end of log.likelihood.MK (Mandelbrot-Koch)


######### putting Likelihood Ratio Test together, June 4, 2015

LRT.MK=function(X1,X2,R=1,w1=1,w2=1,init.mc=0,init.m1=0,init.m2=0,tol=1e-9,show.iter=FALSE) {

Mc=combo.newton.MK(X1,X2,R, w1, w2,init.m=init.mc,tol=tol,show.iter=show.iter)

L0=log.likelihood.MK(X1,Mc,w1)+log.likelihood.MK(X2,R*Mc,w2)

M1=newton.MK(X1,w1,init.m=init.m1,tol=tol,show.iter=show.iter)

M2=newton.MK(X2,w2,init.m=init.m2,tol=tol,show.iter=show.iter)

L1=log.likelihood.MK(X1,M1,w1)+log.likelihood.MK(X2,M2,w2)

chi.sq.stat=-2*(L0-L1)

pval=1-pchisq(chi.sq.stat,1)

return(c(chi.sq.stat,pval))
}




