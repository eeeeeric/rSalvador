######## October 7, 2015: making the B0 distribution work for variation in Nt
### Incorporated into rSalvador 1.4 on December 7, 2015


# ------------------- Wierdl method, incorporated December 9, 2015
### Wierdl estimator, Oct 19, 2015

wierdl.est=function(y,Nt,init.m=0,tol=1e-9,max.iter=25,show.iter=FALSE) {

n=length(y)

lea.est=numeric(n)

for (i in 1:n) {

  lea.est[i]=likely.average(y[i],init.m=init.m, tol=tol,max.iter=max.iter,show.iter=show.iter)/Nt[i]  }

return( median(lea.est) )

} ### --------- end of wierdl.est --------------


# ---------------------- MLE for the benchmark model, when all Nt can be determined ------------

### golden section search for the TMB (theoretical benchmark) model, Aug 2, 2015

golden.benchmark.LD=function(data, Nt, mu.low=1e-12, mu.up=5e-7, tol=1e-9, max.iter=100, show.iter=FALSE) {
      gr = (sqrt(5.0) - 1)/2; 
      a = mu.low
      b = mu.up 
      c = b - gr*(b - a) 
      d = a + gr*(b - a)  
      if (show.iter) message(c(a, b))
      for (i in 1:max.iter){ 
       fc = neg.loglik.LD(data,Nt,c)
       fd = neg.loglik.LD(data,Nt,d)
       if (fc < fd)  { b = d; d = c; c = b - gr*(b - a) }
        else { a = c; c = d; d = a + gr*(b - a)}
      if (show.iter) message( paste("iteration ", toString(i), " ---  [a,b] = [", toString(c(a,b)),"]", sep=''))
       if (abs(b - a) < tol*(abs(c)+abs(d)) ) {return((a + b)/2) }   } ### end of for
}  ### end of golden.LD


### negative of the log likelihood function accouting for every Nt, a theoretical bench mark model, Aug 2, 2015 
neg.loglik.LD=function(y, Nt, b){
     n = length(y)  
     loglikely = 0  
     for(i in 1:n ) {
       m = Nt[i]*b 
       p = tail(prob.LD(m, 1.0, y[i]),1) 
       loglikely = loglikely + log(p) }
      return(-1*loglikely)  }
# --------------- end of the benchmark model ---------------------------


# ------------- beginning of the B0 method ------------------------

# ================= newton.B0, Oct 18, 2015, with CV given as known =======================

newton.B0=function(data, cv=0.1, init.m=0, tol=1e-8, max.iter=30, show.iter=FALSE){

if (init.m>0) {A0=init.m*cv*cv}

else { if (median(data)>0) A0=jones.median.plating(data,e=1)*cv*cv else A0=LD.p0.est(data)*cv*cv }

m0=general.newton.B0(data,k=1/cv/cv, init.A=A0, tol=tol, max.iter=max.iter, show.iter=show.iter)

return(m0/cv/cv) 

}  ### end of newton.B0
 

# ===============  confint.B0, Oct 18, 2015 =====================

confint.B0=function(data,cv=0.1,alpha=0.05,tol=1e-8,init.m=5.1,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

if (init.m>0) {A0=init.m*cv*cv}

else { if (median(data)>0) A0=jones.median.plating(data,e=1)*cv*cv else A0=LD.p0.est(data)*cv*cv }

initLow=init.lower*cv*cv

initUp=init.upper*cv*cv

ci.orig=general.confint.B0(data,k=1/cv/cv,alpha=alpha,init.A=A0,init.lower=initLow,init.upper=initUp,
  max.iter=max.iter, show.iter=show.iter)

return(ci.orig/cv/cv)

}  ######### end of confint.B0


# -------------- Auxiliary functions for B0 distribution ------------


# ===================== newton.B0 October 7, 2015  ===============

general.newton.B0=function(data, k=1, tol=1e-8, init.A=0, max.iter=30, show.iter=FALSE) {

if (init.A>0) {A0=init.A}

else { if (median(data)>0) A0=jones.median.plating(data,e=1)/k else A0=LD.p0.est(data)/k }

if (show.iter) message( paste('iteration 0 yielding ...', toString(A0) ) )

for (i in 1:max.iter) {

score.info=find.score.B0(A0,k,data)

A1=A0+score.info[1]/score.info[2]

if ( abs(A1-A0)/A0 <tol ) {return (A1)}

else {

if (show.iter) message( paste('iteration', toString(i), 'yielding ... ', toString(A1) )) ;

 A0=A1

}  # end of if-else

} # end of for 

return('no convergence')

}



# ===============  confint.B0, Oct 17, 2015 =====================

general.confint.B0=function(data,k=100,alpha=0.05,tol=0.000001,init.A=0.01,init.lower=0,init.upper=0,
  max.iter=30,show.iter=FALSE) {

initLow=init.lower

initUp=init.upper

Ahat=general.newton.B0(data,k=k,tol=tol,init.A=init.A)

if (show.iter) {message( paste('The ML estimate of A ... ',toString(Ahat) ) ) }

score.info.like=find.score.info.loglikely.B0(Ahat,k=k,data=data)

score=score.info.like[1]

info=score.info.like[2]

like=score.info.like[3]

qa=qchisq(1-alpha,1)

h=sqrt(qa/info)

la=like-0.5*qa

if (show.iter) message('Iterating for lower limit ... ')

if (initLow<=0) {A0=Ahat-0.5*h} else {A0=initLow}

for (i in 1:max.iter) {

 like.score=find.loglik.score.B0(A0,k,data)

 like=like.score[1];  score=like.score[2]

 A1=A0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(A1)) ) }

 if ( abs( (A1-A0)/A0 )<tol ) {AL=A1; break} else {A0=A1}

} # end of 1st for loop

if (show.iter) message('Iterating for upper limit ... ')

if (initUp<=0) {A0=Ahat+0.5*h} else {A0=initUp}

for (i in 1:max.iter) {

 like.score=find.loglik.score.B0(A0,k,data)

 like=like.score[1];  score=like.score[2]

 A1=A0 - (like -la)/score

 if (show.iter) {message( paste('iteration ',toString(i), ' yielding', toString(A1)) ) }

 if ( abs( (A1-A0)/A0 )<tol ) {AU=A1; break} else {A0=A1}

} # end of 2nd for loop

return (c(AL,AU))

}  # ============== end of function general.confint.B0 ==============




# ------------- find the score and Fisher information for the B0 model, Oct 7, 2015

find.score.B0=function(A,k,data) {

n=max(data)

p=prob.B0(A,k,n)

deriv=find.B0.deriv(A,k,n)

p1=deriv[[1]]; p2=deriv[[2]]

score=(p1/p)[data+1]

info=( (p1/p)^2-p2/p )[data+1]

return( c( sum(score), sum(info) ) )

}


# ---------- find the first and second derivative, Oct 7, 2015

find.B0.deriv=function(A,k,n) {

h=1/(1:n)/(2:(n+1))

h=append(h,-1,after=0)

h2=seq.convolute(h,h)

p1=k*seq.convolute(prob.B0(A,k+1,n),h)

p2=k*(k+1)*seq.convolute(prob.B0(A,k+2,n),h2)

return( list(p1,p2) )

}


# --------------- log-likelihood and the socre, without Fisher info, Oct 8, 2015

find.loglik.score.B0=function(A,k,data) {

n=max(data)

p=prob.B0(A,k,n)

h=1/(1:n)/(2:(n+1))

h=append(h,-1,after=0)

p1=k*seq.convolute(prob.B0(A,k+1,n),h)

loglik=log( p[data+1] )

score=(p1/p)[data+1]

return( c(sum(loglik),sum(score)) )

}



# --------------- log-likelihood, the socre and the Fisher info: all three, Oct 8, 2015

find.score.info.loglikely.B0=function(A,k,data) {

n=max(data)

p=prob.B0(A,k,n)

h=1/(1:n)/(2:(n+1))

h=append(h,-1,after=0)

h2=seq.convolute(h,h)

p1=k*seq.convolute(prob.B0(A,k+1,n),h)

p2=k*(k+1)*seq.convolute(prob.B0(A,k+2,n),h2)

loglik=log( p[data+1] )

score=(p1/p)[data+1]

info=( (p1/p)^2-p2/p )[data+1]

return( c(sum(score), sum(info), sum(loglik)) )

}
