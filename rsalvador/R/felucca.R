# Nov 9, 2020, tackling fold change for the LD-with plating model

LD.plating.fold.score.info = function(x, y, Nx, Ny, b0, b1, e1 = 0.1, e2 = 0.1) {

  nx = max(x)

  ny = max(y)
   
  mx = Nx * exp(b0)
 
  my = Ny * exp(b0 + b1)

  px = prob.LD.plating(mx, e1, nx)
 
  py = prob.LD.plating(my, e2, ny)

  k1 = e1 / (1 - e1)
   
  k2 = e2 / (1 - e2)

  etaSeq.x = etaSeq(e1, nx)
   
  etaSeq.y = etaSeq(e2, ny)

  px1 = seq.convolute(etaSeq.x, px) * k1  

  px2 = seq.convolute(etaSeq.x, px1) * k1

  py1 = seq.convolute(etaSeq.y, py) * k2  

  py2 = seq.convolute(etaSeq.y, py1) * k2

  F1 = sum((px1 / px)[x + 1])

  F12 = sum(((px1 / px)^2)[x + 1])

  F2 = sum((px2 / px)[x + 1])

  F3 = sum((py1 / py)[y + 1])

  F32 = sum(((py1 / py)^2)[y + 1])

  F4 = sum((py2 / py)[y + 1])

  U2 = Ny * exp(b0 + b1) * F3

  U1 = U2 + Nx * exp(b0) * F1

  J12 = U2 + Ny^2 * exp(2 * (b0 + b1)) * (F4 - F32)

  J22 = J12

  J11 = J12 + Nx * exp(b0) * F1 + Nx^2 * exp(2 * b0) * (F2 - F12)

  return(c(U1, U2, -J11, -J12, -J22))

}


newton.foldchange.LD.plating = function (x, y, Nx, Ny, e1 = 0.1, e2 = 0.1, tol = 1e-08, init.base = 1.5e-8, init.fold = 1.6, 
                                         max.iter = 30, no.log=TRUE, show.iter = FALSE) {   

  LARGE_M = 500; SMALL_M = 1e-9  # to avoid over/under flow, Nov 8, 2020

  beta0=matrix(c(log(init.base), log(init.fold)), ncol=2)

  if (show.iter) 

    message(paste("iteration 0 yielding ... ", toString(beta0)))

  for (i in 1:max.iter) {

  if ((Nx*exp(beta0[1]) > LARGE_M) || (Ny*exp(beta0[1] + beta0[2]) > LARGE_M) || 
      (Nx*exp(beta0[1]) < SMALL_M) || (Ny*exp(beta0[1] + beta0[2]) < SMALL_M))  {

     message('A very large/small value of m has been encountered. New starting value may help.'); return(NA)}  # Nov 8, 202

  srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, beta0[1], beta0[2], e1 = e1, e2 = e2)

  score = matrix(c(srinfo[1], srinfo[2]), ncol=2) 

  fisher = matrix(c(srinfo[3], srinfo[4], srinfo[4], srinfo[5]), ncol=2, byrow=TRUE)   

  beta1 = beta0 + crossprod(t(score), solve(fisher))

    if (sqrt(sum((beta0 - beta1)^2 / sum(beta0^2))) < tol) {
      
      if (no.log) beta1 = exp(beta1)
      
      return(as.vector(beta1)) }

    else {
            if (show.iter) message(paste("iteration", toString(i), "yielding ... ", toString(beta1)))

            beta0 = beta1
        }
    }

  cat('after ', toString(max.iter), ' iterations, no convergence achieved.\n')

  return(NA)

}  # end of newton.foldchange.LD.plating



### ========================== profile CI for b1, mutation rate fold change, Nov 3, 2020 ===========

# ----------- define starting value for constructing CI for b1, Nov 3, 2010 ----------

fold.LD.plating.venzon.b1 = function(x, y, Nx, Ny, b0.hat, b1.hat, alpha=0.05, e1 = 0.1, e2 = 0.1) {

  srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, b0.hat, b1.hat, e1 = e1, e2 = e2)

  q = qchisq(1 - alpha, df = 1)

  J11 = srinfo[3]; J12 = srinfo[4]; J22 = srinfo[5]

  h = J22 - J12^2 / J11

  h = 0.5 * sqrt(q / h)

  ven.low.b0 = b0.hat + h * J12 / J11; ven.low.b1 = b1.hat - h

  ven.up.b0 = b0.hat - h * J12 / J11; ven.up.b1 = b1.hat + h

return( c(ven.low.b0, ven.low.b1, ven.up.b0, ven.up.b1) ) 
} 


# ----------- the log likelihood needed for calculating init values, Nov 3, 2020 ---------------

fold.loglik.LD.plating = function(x, y, Nx, Ny, b0, b1, e1 = 0.1, e2 = 0.1) {

  mx = Nx*exp(b0)

  my = Ny * exp(b0 + b1)

  probx = prob.LD.plating(m = mx, e = e1, n = max(x))

  proby = prob.LD.plating(m = my, e = e2, n = max(y))

  loglik = sum(log(probx[x + 1])) + sum(log(proby[y + 1]))

return(loglik)

}


# --------------- computing profile CI for b1, the logarithm of fold change, Nov 3, 2020 ---------------

confint.foldchange.LD.plating = function(x, y, Nx, Ny, e1 = 0.1, e2 = 0.1, alpha = 0.05, init.base = 1e-8,
                     init.fold = 1.6, init.low.base = -9, init.up.base = -9, init.low.fold = -9,
                     init.up.fold = -9, max.iter = 30, tol = 1e-9, show.iter = FALSE) {

  init.b00 = log(init.base)

  init.b10 = log(init.fold)

  if (show.iter) cat('\nComputing MLE ...\n')

  mle = newton.foldchange.LD.plating(x, y, Nx, Ny, e1 = e1, e2 = e2, init.base = init.base, init.fold = init.fold, 
                                 tol = tol, max.iter = max.iter, no.log = FALSE, show.iter = show.iter)

  if (is.na(mle)[1]) return(NA) # Nov 8, 2020

  if (show.iter) cat('MLE of b0 and b1 ...',  mle, ' \n')

  b0.hat = mle[1]; b1.hat = mle[2]

  if ((init.low.base < 0) || (init.low.fold < 0) || (init.up.base < 0) || (init.up.base < 0))

    {all.inits = fold.LD.plating.venzon.b1(x, y, Nx, Ny, b0.hat, b1.hat, e1 = e1, e2 = e2, alpha = alpha) }

  if (init.low.base < 0) {b00 = all.inits[1]} else b00 = log(init.low.base)

  if (init.low.fold < 0) {b10 = all.inits[2]} else b10 = log(init.low.fold)


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

  srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, b00, b10, e1 = e1, e2 = e2)

  zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

  like = fold.loglik.LD.plating(x, y, Nx, Ny, b00, b10, e1 = e1, e2 = e2)

  la = fold.loglik.LD.plating(x, y, Nx, Ny, b0.hat, b1.hat, e1 = e1, e2 = e2) - 0.5 * qchisq(1 - alpha, df = 1)

  zheng37vet=matrix(c(like - la, srinfo[1]), ncol = 1)

  beta0 = matrix(c(b00,b10), ncol = 1)

  if(show.iter) cat('\nComputing the lower limit ... \n')

  if (show.iter) cat("Iteration ", 0, " ... ", beta0[2], " \n")

  for(i in 1:max.iter) {

    if (i == max.iter) {message('Maximum number of iterations reached.'); return(NA)} # Nov 10, 2020

    beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

    if (show.iter) cat("Iteration ", i, " ... ", beta1[2], " \n")

  ### updating .........

    b01 = beta1[1]; b11 = beta1[2]

    b00 = beta0[1]; b10 = beta0[2]

    if (sqrt(sum((beta0 - beta1)^2 / sum(beta0^2))) < tol) { b1.lower = b11; break() }

    srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, b01, b11, e1 = e1, e2 = e2)

    zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

    like=fold.loglik.LD.plating(x, y, Nx, Ny, b01, b11, e1 = e1, e2 = e2)

    zheng37vet=matrix(c(like - la, srinfo[1]), ncol = 1)

    beta0 = beta1

    } ## end of for lower limit

# --------------- repeating for the upper limit --------------------


  if (init.up.base<0) {b00 = all.inits[3]} else b00 = log(init.up.base)

  if (init.up.fold<0) {b10 = all.inits[4]} else b10 = log(init.up.fold)

### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

  srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, b00, b10, e1 = e1, e2 = e2)

  zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

  like = fold.loglik.LD.plating(x, y, Nx, Ny, b00, b10, e1 = e1, e2 = e2)

  la = fold.loglik.LD.plating(x, y, Nx, Ny, b0.hat, b1.hat, e1 = e1, e2 = e2) - 0.5 * qchisq(1 - alpha, df = 1)

  zheng37vet = matrix(c(like - la, srinfo[1]), ncol = 1)

  beta0 = matrix(c(b00, b10), ncol = 1)

  if (show.iter) cat('\nComputing the upper limit ... \n')

  if (show.iter) cat("Iteration ", 0, " ... ", beta0[2], " \n")

  for(i in 1:max.iter) {

    if (i == max.iter) {message('Maximum number of iterations reached.'); return(NA)} # Nov 10, 2020

    beta1 = beta0 + solve(zheng37mat) %*% zheng37vet

    if (show.iter) cat("Iteration ", i, " ... ", beta1[2], " \n")

### updating .........

    b01 = beta1[1]; b11 = beta1[2]

    b00 = beta0[1]; b10 = beta0[2]

    if (sqrt(sum((beta0 - beta1)^2 / sum(beta0^2))) < tol) { b1.upper = b11; break() }

    srinfo = LD.plating.fold.score.info(x, y, Nx, Ny, b01, b11, e1 = e1, e2 = e2)

    zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

    like = fold.loglik.LD.plating(x, y, Nx, Ny, b01, b11, e1 = e1, e2 = e2)

    zheng37vet = matrix(c(like - la, srinfo[1]), ncol = 1)

    beta0 = beta1

    } ## end of for upper limit

  if (show.iter) cat('\n')

  return(exp(c(b1.lower, mle[2], b1.upper)))

}  ########## end of confint.foldchange.LD, Nov 4, 2020



########### Bayes ###########
# -------- Bayesian approach to fold change, based on existing mcmc.LD.plating.test, Nov 14, 2020 -------------

bayes.foldchange.LD.plating=function(y1, y2, N1, N2, e1 = 0.1, e2 = 0.1, s0 = 0.4, s1 = 0.5, 
                                     init.base = 1e-8, init.fold = 1.6, v0 = 100, v1 = 100,
                                     iter = 1100, burn = 100, thin = 1, alpha = 0.05, short.out = TRUE, show.simu = FALSE) {

  mu0 = log(init.base)

  mu1 = log(init.fold)

  q1 = alpha/2;  q2=1-q1;

  simu.chain = mcmc.LD.plating.test(y1, y2, N1, N2, e1 = e1, e2 = e2, s0 = s0, s1 = s1, mu0 = mu0, mu1 = mu1,

       v0=v0, v1 = v1, Iter = iter, show.simu = show.simu)

  what.take = seq(burn + 1, iter, by = thin)

  chain1 = simu.chain[[2]][what.take, 1]

  chain2 = simu.chain[[2]][what.take, 2]

  out1 = exp(quantile(chain1, c(q1, 0.5, q2)))

  out2 = exp(quantile(chain2, c(q1, 0.5, q2)))

  accept.rate = simu.chain[[1]]

  if (short.out) out.stuff = out2

  else out.stuff = list(accept.rate, list(out1, out2), list(chain1, chain2))

  return(out.stuff)

}


# ------------------ Bootstrap approach, Nov 7, 2020 -------------

boot.LD.plating.one = function(y1, y2, N1, N2, e1 = 0.1, e2 = 0.1, init.m1 = 0, init.m2 = 0) {

  n1 = length(y1)

  n2 = length(y2)

  z1 = sample(y1, n1, replace = TRUE)

  z2 = sample(y2, n2, replace = TRUE)

  mu1 = newton.LD.plating(z1, e = e1, init.m = init.m1) / N1

  mu2 = newton.LD.plating(z2, e = e2, init.m = init.m2) / N2

  return(mu2 / mu1)

}


boot.foldchange.LD.plating = function(y1, y2, N1, N2, e1 = 0.1, e2 = 0.1, alpha = 0.05, init.m1 = 0, init.m2 = 0, 
   na.ok = FALSE, nboot = 10) {
 
  q1 = alpha / 2; q2 = 1 - q1;

  rates = rep(-1, nboot)

  for (i in 1:nboot) rates[i] = boot.LD.plating.one(y1, y2, N1, N2, e1 = e1, e2 = e2, init.m1 = init.m1, init.m2 = init.m2)

  n.fail=sum(sapply(rates, is.na))
 
  if (n.fail == 0) return(quantile(rates, c(q1, 0.5, q2)))

    else if (na.ok) {message(paste(toString(n.fail), 'bootstrap samples out of', toString(nboot), 'did not allow
      convergence when calculating mutation rates. New starting values may be helpful.' ));  
        return(quantile(rates, c(q1, 0.5, q2), na.rm = TRUE))}

    else return (NA)
}


