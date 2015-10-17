# rSalvador: An R tool for the Luria-Delbruck fluctuation assay
# Qi Zheng, Department of epidemiology and Biostatistics
# Texas A&M School of Public Health
# Version 1.0: April 20, 2014
# Version 1.1: April 19, 2015
# Version 1.2: June 2, 2015
# Version 1.3: June 28, 2015
# For non-essential functions file, June 5, 2015

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

### --------------------- likelihood ratio test for LD model ------------------
# April 18, 2015 -- to solve Dr. Feher's comparison problem
# May 6, 2015 -- combine both functions for reviewers convenience

### finding the score and observed info for the model under null hypothesis

combo.LD.score.info=function(m,X1,X2,R, phi1=1,phi2=1) {

us1=LD.score.info(m,phi1,X1)

us2=LD.score.info(R*m,phi2,X2)

u1=us1[1]

j1=us1[2]

u2=us2[1]

j2=us2[2]

u=u1+R*u2

j=j1+R*R*j2

return(c(u,j))

}

### find Mc, the combined MLE, April 19, 2015


# ================== combind newtonLD ===============

combo.newton.LD=function(X1,X2,R, phi1=1.0, phi2=1.0,
           tol=0.00000001, init.m=0, max.iter=30, show.iter=FALSE) {

Xc=c(X1,X2)

if (init.m>0) {m0=init.m}

   else if (init.m==0) 

       {if (min(Xc)==0) m0=LD.p0.est(Xc) else m0=jones.median.est(Xc) }  ## option A

   else 

       {if (min(X1)==0) m0=LD.p0.est(X1) else m0=jones.median.est(X1) }  ## option B

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

for (i in 1:max.iter) {

score.info=combo.LD.score.info(m0,X1,X2,R,phi1,phi2)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}


######### putting all together, April 20, 2015

LRT.LD=function(X1,X2,R=1,phi1=1,phi2=1,init.mc=0,init.m1=0,init.m2=0,tol=1e-9,show.iter=FALSE) {

Mc=combo.newton.LD(X1,X2,R, phi1, phi2,init.m=init.mc,tol=tol,show.iter=show.iter)

L0=log.likelihood.LD(X1,Mc,phi1)+log.likelihood.LD(X2,R*Mc,phi2)

M1=newton.LD(X1,phi1,init.m=init.m1,tol=tol,show.iter=show.iter)

M2=newton.LD(X2,phi2,init.m=init.m2,tol=tol,show.iter=show.iter)

L1=log.likelihood.LD(X1,M1,phi1)+log.likelihood.LD(X2,M2,phi2)

chi.sq.stat=-2*(L0-L1)

pval=1-pchisq(chi.sq.stat,1)

return(c(chi.sq.stat,pval))
}

### ----------- end of LRT.LD -----------------------





###############  LTR.LD.plating, likelihood ratio test  ####################

# May 6, 2015, small modification June 29, 2015 

### finding the score and observed info for the model under null hypothesis

combo.LD.score.info.plating=function(m,X1,X2,R, e1=0.1,e2=0.1) {

eta1=etaSeq(e1,max(X1))

eta2=etaSeq(e2,max(X2))

us1=LD.score.info.plating(m,e1,eta1,X1)

us2=LD.score.info.plating(R*m,e2,eta2,X2)

u1=us1[1]

j1=us1[2]

u2=us2[1]

j2=us2[2]

u=u1+R*u2

j=j1+R*R*j2

return(c(u,j))

}

### find Mc, the combined MLE with plating, April 22, 2015

combo.newton.LD.plating=function(X1,X2,R,e1,e2, tol=0.00000001, init.m=0, max.iter=30,
    show.iter=FALSE) {

Xc=c(X1,X2)

n1=length(X1); n2=length(X2); ec=(n1*e1+n2*e2)/(n1+n2)

if (init.m>0) {m0=init.m}

    else if (init.m==0) 
   
         {if (median(Xc)==0) m0=p0.plating(Xc,ec) else m0=jones.median.plating(Xc,ec) } #option A

    else 

        {if (min(X1)==0) m0=p0.plating(X1,e1) else m0=jones.median.plating(X1,e1) }  #option B

if (show.iter) message( paste('iteration 0 yielding ...', toString(m0) ) )

eta1=etaSeq(e1,max(X1))

eta2=etaSeq(e2,max(X2))

for (i in 1:max.iter) {

score.info=combo.LD.score.info.plating(m0,X1,X2,R,e1,e2)

m1=m0+score.info[1]/score.info[2]

if ( abs(m1-m0)/m0 <tol ) {return (m1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(m1) )) ;

 m0=m1

}  # end of if-else

} # end of for 

return('no convergence')

}



######### putting all together -- comparison with plating, April 22, 2015

LRT.LD.plating=function(X1,X2,R=1,e1=0.1,e2=0.1,init.mc=0,init.m1=0,init.m2=0,
     tol=1e-9,show.iter=FALSE) {

Mc=combo.newton.LD.plating(X1,X2,R,e1,e2,init.m=init.mc,tol=tol,show.iter=show.iter)

L0=log.likelihood.LD.plating(X1,Mc,e1)+log.likelihood.LD.plating(X2,R*Mc,e2)

M1=newton.LD.plating(X1,e1,init.m=init.m1,tol=tol,show.iter=show.iter)

M2=newton.LD.plating(X2,e2,init.m=init.m2,tol=tol,show.iter=show.iter)

L1=log.likelihood.LD.plating(X1,M1,e1)+log.likelihood.LD.plating(X2,M2,e2)

chi.sq.stat=-2*(L0-L1)

pval=1-pchisq(chi.sq.stat,1)

return(c(chi.sq.stat,pval))
}


#####################################################
### ----------- MCMC --------------------
##### MCMC for the MK model, June 28, 2015

mcmc.MK=function(y,Nt,w=1,Iter=100,init.mu=2e-8,b0=-15,v0=30,s0=0.8,show.simu=FALSE) {

simu.b=matrix(nrow=Iter,ncol=1)

current.b=log(init.mu)

tot.acc=0

for (i in 1:Iter) {

prop.b=rnorm(1,current.b,sqrt(s0))

h1=log.likelihood.MK(y,Nt*exp(prop.b),w)-(prop.b-b0)^2/(2*v0)+prop.b  ### added b 7-1-2015

h2=log.likelihood.MK(y,Nt*exp(current.b),w)-(current.b-b0)^2/(2*v0)+current.b  ### 7-1-2015

h=h1-h2

u=log(runif(1))

if (u<h) {current.b=prop.b; tot.acc=tot.acc+1} 

simu.b[i]=current.b

if (show.simu) {print(current.b)}

}  ### end of loop

return(list(tot.acc/Iter,simu.b))
}  ### end of mcmc.MK

### ---------- mutation rate comparison by MCMC ---------------
# computing log-likelihood to be used in the posterior distribution, June 28, 2015

posterior.MK=function(x,y,Nx,Ny,w1=1,w2=1,b0,b1) {

m1=Nx*exp(b0)

m2=Ny*exp(b0+b1)

h1=log.likelihood.MK(x,m1,w=w1)

h2=log.likelihood.MK(y,m2,w=w2)

h=h1+h2+2*b0+b1   ### adding b0+b1 as log of Jacobian, July 1, 2015 and July 2, 2015

return(h)  } ##### end of posterior.MK



### ----------- mcmc.MK.test, June 28, 2015

mcmc.MK.test=function(x,y,N1,N2,w1=1,w2=1,Iter=100,init.mu0=-20, init.mu1=0.05, mu0=-20.7,
                              mu1=2.7,v0=30,v1=30,s0=0.5,s1=0.5,show.simu=FALSE){

beta=matrix(nrow=Iter,ncol=2)

current.beta=c(init.mu0,init.mu1)

tot.acc=c(0,0)

prop.s=c(s0,s1)   ### tuning parameters

for (t in 1:Iter) {

prop.beta=current.beta

prop.beta[1]=rnorm(1,current.beta[1],prop.s[1])

h1=posterior.MK(x,y,N1,N2,w1,w2,prop.beta[1],prop.beta[2])-
   posterior.MK(x,y,N1,N2,w1,w2,current.beta[1],current.beta[2])

h2=  -(prop.beta[1]-mu0)^2/(2*v0)-(prop.beta[2]-mu1)^2/(2*v1)
  +(current.beta[1]-mu0)^2/(2*v0)+(current.beta[2]-mu1)^2/(2*v1)

loga=h1+h2

u=log(runif(1))

if (u<loga) {current.beta=prop.beta;tot.acc[1]=tot.acc[1] +1}

prop.beta=current.beta

prop.beta[2]=rnorm(1,current.beta[2],prop.s[2])

h1=posterior.MK(x,y,N1,N2,w1,w2,prop.beta[1],prop.beta[2])-
   posterior.MK(x,y,N1,N2,w1,w2,current.beta[1],current.beta[2])

h2=  -(prop.beta[1]-mu0)^2/(2*v0)-(prop.beta[2]-mu1)^2/(2*v1)
  +(current.beta[1]-mu0)^2/(2*v0)+(current.beta[2]-mu1)^2/(2*v1)

loga=h1+h2

u=log(runif(1))

if (u<loga) {current.beta=prop.beta; tot.acc[2]=tot.acc[2]+1}

beta[t,]=current.beta

if (show.simu) {print(current.beta)}

}  ### end of simulation loop

return( list(tot.acc/Iter,beta) )

}  ### end of mcmc.MK.test


##### MCMC for the LD model with partial plating, June 29, 2015

mcmc.LD.plating=function(y,Nt,e=0.1,Iter=100,init.mu=2e-8,b0=-15,v0=30,s0=0.8,show.simu=FALSE) {

simu.b=matrix(nrow=Iter,ncol=1)

current.b=log(init.mu)

tot.acc=0

for (i in 1:Iter) {

prop.b=rnorm(1,current.b,sqrt(s0))

h1=log.likelihood.LD.plating(y,Nt*exp(prop.b),e)-(prop.b-b0)^2/(2*v0)+prop.b  ## 7-1-2015

h2=log.likelihood.LD.plating(y,Nt*exp(current.b),e)-(current.b-b0)^2/(2*v0)+current.b  ## 7-1-2015

h=h1-h2

u=log(runif(1))

if (u<h) {current.b=prop.b; tot.acc=tot.acc+1} 

simu.b[i]=current.b

if (show.simu) {print(current.b)}

}  ### end of loop

return(list(tot.acc/Iter,simu.b))
}  ### end of mcmc.LD.plating: June 29, 2015

### ---------- mcmc comparison under LD model with partial plating, June 29, 2015

# computing log-likelihood to be used in the posterior distribution, June 29, 2015

posterior.LD.plating=function(x,y,Nx,Ny,e1=0.1,e2=0.1,b0,b1) {

m1=Nx*exp(b0)

m2=Ny*exp(b0+b1)

h1=log.likelihood.LD.plating(x,m1,e=e1)

h2=log.likelihood.LD.plating(y,m2,e=e2)

h=h1+h2+2*b0+b1  ### July 1, 2015: add b0+b1 as the Jacobian and July 2, 2015

return(h)  } ##### end of posterior.LD.plating



### ----------- mcmc.LD.plating.test, June 29, 2015

mcmc.LD.plating.test=function(x,y,N1,N2,e1=0.1,e2=0.1,Iter=100,init.mu0=-20,
    init.mu1=0.05, mu0=-20.7, mu1=2.7,v0=30,v1=30,s0=0.5,s1=0.5,show.simu=FALSE){

beta=matrix(nrow=Iter,ncol=2)

current.beta=c(init.mu0,init.mu1)

tot.acc=c(0,0)

prop.s=c(s0,s1)   ### tuning parameters

for (t in 1:Iter) {

prop.beta=current.beta

prop.beta[1]=rnorm(1,current.beta[1],prop.s[1])

h1=posterior.LD.plating(x,y,N1,N2,e1,e2,prop.beta[1],prop.beta[2])-
   posterior.LD.plating(x,y,N1,N2,e1,e2,current.beta[1],current.beta[2])

h2=  -(prop.beta[1]-mu0)^2/(2*v0)-(prop.beta[2]-mu1)^2/(2*v1)
  +(current.beta[1]-mu0)^2/(2*v0)+(current.beta[2]-mu1)^2/(2*v1)

loga=h1+h2

u=log(runif(1))

if (u<loga) {current.beta=prop.beta;tot.acc[1]=tot.acc[1] +1}

prop.beta=current.beta

prop.beta[2]=rnorm(1,current.beta[2],prop.s[2])

h1=posterior.LD.plating(x,y,N1,N2,e1,e2,prop.beta[1],prop.beta[2])-
   posterior.LD.plating(x,y,N1,N2,e1,e2,current.beta[1],current.beta[2])

h2=  -(prop.beta[1]-mu0)^2/(2*v0)-(prop.beta[2]-mu1)^2/(2*v1)
  +(current.beta[1]-mu0)^2/(2*v0)+(current.beta[2]-mu1)^2/(2*v1)

loga=h1+h2

u=log(runif(1))

if (u<loga) {current.beta=prop.beta; tot.acc[2]=tot.acc[2]+1}

beta[t,]=current.beta

if (show.simu) {print(current.beta)}

}  ### end of simulation loop

return( list(tot.acc/Iter,beta) )

}  ### end of mcmc.LD.plating.test



####### August 2, 2015: the B0 Bartlett distribution

prob.B0=function(A,k,n=5) {

result.p=rep(0,n+1)

### if (! is.loaded('bart0c')) dyn.load('bart0.so')

z=.C('pmfBartzero_R_wrapper', A=as.double(A), k=as.double(k), n=as.integer(n), prob=as.double(result.p))

return(z$prob)

}


# ------------ accouting for variation in Nt under LD model, added Aug 3, 2015 ----------
### Revised, October 4, 2015 -- deal with m directly, not via log(m)

### negative of log likelihood function of the Bartlett B0 model, Aug 2, 2015
neg.loglik.B0=function(data, cv, m0) {
     nn = max(data);
     p = prob.B0(cv*cv*m0, 1/cv/cv, nn) 
     lik = sum(log(p[data + 1]))
     return(-1*lik)  }

### golden section search for the Bartlett B0 model, cv is known, Aug 2, 2015
golden.LD.B0=function(data, cv=0.1, m.low=0.1, m.up=30, tol=1e-8, max.iter=60, show.iter=FALSE) {
      gr = (sqrt(5.0) - 1)/2; 
      a = m.low
      b = m.up 
      c = b - gr*(b - a) 
      d = a + gr*(b - a)  
      if (show.iter)
         message( paste("iteration ", toString(0), " ---  (a,b) = (", toString(c(a,b)),")"))
      for (i in 1:max.iter){ 
       fc = neg.loglik.B0(data,cv,c)
       fd = neg.loglik.B0(data,cv,d)
       if (fc < fd)  { b = d; d = c; c = b - gr*(b - a) }
        else { a = c; c = d; d = a + gr*(b - a)}
      if (show.iter)
         message( paste("iteration ", toString(i), " ---  (a,b) = (", toString(c(a,b)),")"))
       if ( abs(a - b) < tol*(abs(c)+abs(d)) ) {return((a + b)/2) } } ### end of for
}  ### end of golden.LD.B0, modified Oct 5, 2015

