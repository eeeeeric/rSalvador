### starting date: Jan 13, 2016, joint inference using the MK model, for rSalvador version 1.5

# -------------- first and second order derivatives, Jan 13, 2016 -------------------

MK.partials=function(m,w,n) {

r=1/w

B=( 2:n -1)/(2:n +r)

B=Reduce(prod, B, init=1,accumulate=TRUE)/(1+r)

eta=1/(1:n +r)

eta=Reduce(sum,eta[-1],init=eta[1],accumulate=TRUE)

h10=append(r*B,-1,after=0)

h01=append(B*(1-r*eta),0,after=0)

p=prob.MK(m=m,w=w,n=n)

p10=seq.convolute(p,h10)

p01=m*seq.convolute(p,h01)

h02=(1:n)/(1:n + r)^2

h02=Reduce(sum,h02[-1],init=h02[1],accumulate=TRUE)

h02=B*(r*eta*eta - eta -h02)

h02=append(h02, 0, after=0)

p02=m*(seq.convolute(h02,p)+seq.convolute(h01,p01))

p20=seq.convolute(h10,p10)

p11=p01/m + m*seq.convolute(h01,p10)

# in the Zheng (2005) paper, partial derivatives were given in terms of r=1/w. Now change r to w ...

q01=-p01/w^2

q11=-p11/w^2

q02=2*p01/w^3 +p02/w^4

return( list(p10, q01, p20, q11, q02, p) )

}

# -------------- score and observed information for the MK model, Jan 14, 2016 ---------

MK.score.info.joint=function(y, m, w) {

n=max(y)

jumble=MK.partials(m, w, n)

p10=jumble[[1]]; p01=jumble[[2]]; p20=jumble[[3]]

p11=jumble[[4]]; p02=jumble[[5]]; p=jumble[[6]]

U1=sum( (p10/p)[y+1] )

U2=sum( (p01/p)[y+1] )

J11=sum( ( (p10/p)^2 -p20/p )[y+1] )

J12=sum( ( p10*p01/p^2 -p11/p)[y+1] )

J22=sum( ( (p01/p)^2 -p02/p)[y+1] )

return(list(U1, U2, J11, J12, J22))

}


# -------------- Newton-Raphson algorithm, Jan 14, 2016 --------------

newton.joint.MK=function (data, tol = 1e-08, init.m = -1, init.w=-1, max.iter = 30, 
    show.iter = FALSE) 
{   

if (init.m < 0) { m0=newton.LD(data) } else {m0=init.m}

if (init.w < 0) { w0=1.0 } else {w0=init.w}

beta0=matrix( c(m0, w0), ncol=2)

if (show.iter) 
        message(paste("iteration 0 yielding ...", toString(beta0)))

for (i in 1:max.iter) {

srinfo=MK.score.info.joint(data, beta0[1], beta0[2] )

score=matrix( c(srinfo[[1]], srinfo[[2]]), ncol=2) 

fisher=matrix(c(srinfo[[3]],srinfo[[4]],srinfo[[4]],srinfo[[5]]), ncol=2, byrow=TRUE)   

beta1 = beta0 + crossprod(t(score), solve(fisher) )

        if ( sqrt(sum((beta0-beta1)^2))  < tol) {
            return(as.vector(beta1))
        }
        else {
            if (show.iter) 
                message(paste("iteration", toString(i), "yielding ... ", 
                  toString(beta1)))
            beta0=beta1
        }
    }
    cat("after ", toString(max.iter), " iterations, no convergence achieved\n.")
}


# --- Jan 15, 2016; profile likelihood-based confidence intervals ---

# ---------- log-likelihood function ------------

loglikMK=function(y,m,w) {

prob=prob.MK(m,w,n=max(y))

loglik=sum(log(prob[y+1]))

return(loglik) }


# ----------- define starting value, Jan 15, 2016 --------------

venzon.m=function(data, m.hat, w.hat, alpha=0.05) {

srinfo=MK.score.info.joint(data, m.hat, w.hat)

q=qchisq(1-alpha, df=1)

J11=srinfo[[3]]; J12=srinfo[[4]]; J22=srinfo[[5]]

h=J11-J12^2/J22

h=0.5*sqrt(q/h)

ven.low.m=m.hat-h; ven.low.w=w.hat+h*J12/J22 

ven.up.m=m.hat+h; ven.up.w=w.hat-h*J12/J22

return( c(ven.low.m, ven.low.w, ven.up.m, ven.up.w) )

} 


# --------------- computing profile CI for m, Jan 15, 2016 ---------------

confint.profile.m=function(data, alpha=0.05, init.low.m=-1, init.low.w=-1, init.up.m=-1,
              init.up.w=-1, init.m=-1, init.w=-1,  max.iter=30, tol=1e-6, show.iter=FALSE) {

if (init.m<0) {init.m0=-1} else {init.m0=init.m}

if (init.w<0) {init.w0=-1} else {init.w0=init.w}

if (show.iter) cat('\nComputing MLE ...\n')

mle=newton.joint.MK(data, init.m=init.m0, init.w=init.w0, tol=tol, max.iter=max.iter)

if (show.iter) cat('MLE of m and w ...',  mle, ' \n')

m.hat=mle[1]; w.hat=mle[2]

if ( (init.low.m<0)||(init.low.w<0)||(init.up.m<0)||(init.up.w<0) )

    {all.inits=venzon.m(data, m.hat, w.hat, alpha=alpha) }

if (init.low.m<0) {m0=all.inits[1]} else m0=init.low.m

if (init.low.w<0) {w0=all.inits[2]}  else w0=init.low.w


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

srinfo=MK.score.info.joint(data,m0,w0)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[4]],srinfo[[5]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m0,w0)

la=loglikMK(data,m.hat,w.hat)-0.5*qchisq(1-alpha,df=1)

zheng37vet=matrix(c(like-la,srinfo[[2]]), ncol=1)

beta0=matrix(c(m0,w0),ncol=1)

if(show.iter) cat('\nComputing the lower limit ... \n')

if (show.iter) cat("Iteration ", 0, " ... ", beta0[1], " \n")

for(i in 1:max.iter) {

beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

if (show.iter) cat("Iteration ", i, " ... ", beta1[1], " \n")

### updating .........

m1=beta1[1]; w1=beta1[2]
m0=beta0[1]; w0=beta0[2]

if ( abs(m1-m0)/m0<tol ) {m.lower=m1; break() }

srinfo=MK.score.info.joint(data,m1,w1)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[4]],srinfo[[5]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m1,w1)

zheng37vet=matrix(c(like-la, srinfo[[2]]), ncol=1)

beta0=beta1

} ## end of for lower limit

# --------------- repeating for the upper limit --------------------


if (init.up.m<0) {m0=all.inits[3]} else m0=init.up.m

if (init.up.w<0) {w0=all.inits[4]}  else w0=init.up.w


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

srinfo=MK.score.info.joint(data,m0,w0)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[4]],srinfo[[5]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m0,w0)

la=loglikMK(data,m.hat,w.hat)-0.5*qchisq(1-alpha,df=1)

zheng37vet=matrix(c(like-la,srinfo[[2]]), ncol=1)

beta0=matrix(c(m0,w0),ncol=1)

if(show.iter) cat('\nComputing the upper limit ... \n')

if (show.iter) cat("Iteration ", 0, " ... ", beta0[1], " \n")

for(i in 1:max.iter) {

beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

if (show.iter) cat("Iteration ", i, " ... ", beta1[1], " \n")

### updating .........

m1=beta1[1]; w1=beta1[2]
m0=beta0[1]; w0=beta0[2]

if ( abs(m1-m0)/m0<tol ) {m.upper=m1; break() }

srinfo=MK.score.info.joint(data,m1,w1)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[4]],srinfo[[5]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m1,w1)

zheng37vet=matrix(c(like-la, srinfo[[2]]), ncol=1)

beta0=beta1

} ## end of for upper limit

if (show.iter) cat('\n')
return( c(m.lower, m.upper) )

}


### ========================== profile CI for w, relative fitness, Jan 16, 2016 ===========

# ----------- define starting value for constructing CI for w, Jan 16, 2016 --------------

venzon.w=function(data, m.hat, w.hat, alpha=0.05) {

srinfo=MK.score.info.joint(data, m.hat, w.hat)

q=qchisq(1-alpha, df=1)

J11=srinfo[[3]]; J12=srinfo[[4]]; J22=srinfo[[5]]

h=J22-J12^2/J11

h=0.5*sqrt(q/h)

ven.low.m=m.hat+h*J12/J11; ven.low.w=w.hat-h

ven.up.m=m.hat-h*J12/J11; ven.up.w=w.hat+h

return( c(ven.low.m, ven.low.w, ven.up.m, ven.up.w) ) 
} 


# --------------- computing profile CI for w, Jan 16, 2016 ---------------

confint.profile.w=function(data, alpha=0.05, init.low.m=-1, init.low.w=-1, init.up.m=-1,
              init.up.w=-1, init.m=-1, init.w=-1, max.iter=30, tol=1e-6, show.iter=FALSE) {

if (init.m<0) {init.m0=-1} else {init.m0=init.m}

if (init.w<0) {init.w0=-1} else {init.w0=init.w}

if (show.iter) cat('\nComputing MLE ...\n')

mle=newton.joint.MK(data, init.m=init.m0, init.w=init.w0, tol=tol, max.iter=max.iter)

if (show.iter) cat('MLE of m and w ...',  mle, ' \n')

m.hat=mle[1]; w.hat=mle[2]

if ( (init.low.m<0)||(init.low.w<0)||(init.up.m<0)||(init.up.w<0) )

    {all.inits=venzon.w(data, m.hat, w.hat, alpha=alpha) }

if (init.low.m<0) {m0=all.inits[1]} else m0=init.low.m

if (init.low.w<0) {w0=all.inits[2]}  else w0=init.low.w


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

srinfo=MK.score.info.joint(data,m0,w0)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[3]],srinfo[[4]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m0,w0)

la=loglikMK(data,m.hat,w.hat)-0.5*qchisq(1-alpha,df=1)

zheng37vet=matrix(c(like-la,srinfo[[1]]), ncol=1)

beta0=matrix(c(m0,w0),ncol=1)

if(show.iter) cat('\nComputing the lower limit ... \n')

if (show.iter) cat("Iteration ", 0, " ... ", beta0[2], " \n")

for(i in 1:max.iter) {

beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

if (show.iter) cat("Iteration ", i, " ... ", beta1[2], " \n")

### updating .........

m1=beta1[1]; w1=beta1[2]
m0=beta0[1]; w0=beta0[2]

if ( abs(w1-w0)/w0<tol ) {w.lower=w1; break() }

srinfo=MK.score.info.joint(data,m1,w1)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[3]],srinfo[[4]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m1,w1)

zheng37vet=matrix(c(like-la, srinfo[[1]]), ncol=1)

beta0=beta1

} ## end of for lower limit

# --------------- repeating for the upper limit --------------------


if (init.up.m<0) {m0=all.inits[3]} else m0=init.up.m

if (init.up.w<0) {w0=all.inits[4]}  else w0=init.up.w


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

srinfo=MK.score.info.joint(data,m0,w0)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[3]],srinfo[[4]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m0,w0)

la=loglikMK(data,m.hat,w.hat)-0.5*qchisq(1-alpha,df=1)

zheng37vet=matrix(c(like-la,srinfo[[1]]), ncol=1)

beta0=matrix(c(m0,w0),ncol=1)

if(show.iter) cat('\nComputing the upper limit ... \n')

if (show.iter) cat("Iteration ", 0, " ... ", beta0[2], " \n")

for(i in 1:max.iter) {

beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

if (show.iter) cat("Iteration ", i, " ... ", beta1[2], " \n")

### updating .........

m1=beta1[1]; w1=beta1[2]
m0=beta0[1]; w0=beta0[2]

if ( abs(m1-m0)/m0<tol ) {w.upper=w1; break() }

srinfo=MK.score.info.joint(data,m1,w1)

zheng37mat=matrix(c(-srinfo[[1]],-srinfo[[2]],srinfo[[3]],srinfo[[4]]), ncol=2, byrow=TRUE)

like=loglikMK(data,m1,w1)

zheng37vet=matrix(c(like-la, srinfo[[1]]), ncol=1)

beta0=beta1

} ## end of for upper limit

if (show.iter) cat('\n')
return( c(w.lower, w.upper) )

}  ########## end of confint.profile.w (Jan 16, 2016)

