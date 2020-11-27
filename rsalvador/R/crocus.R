# Nov 12, 2020, tackling fold change for the MK model


######## May 15, 2015: fitness in Luria-Delbruck

# ----------------------- end of old stuff --------------------



MK.fold.score.info = function(x, y, Nx, Ny, b0, b1, w1 = 0.8, w2 = 0.8) {

  nx = max(x)

  ny = max(y)
   
  mx = Nx * exp(b0)
 
  my = Ny * exp(b0 + b1)

  px = prob.MK(mx, w1, nx)
 
  py = prob.MK(my, w2, ny)

  Bx = betaSeq(w1, nx) / w1
  
  By = betaSeq(w2, ny) / w2

  hx = append(Bx, -1, after = 0)  

  hy = append(By, -1, after = 0)  

  px1 = seq.convolute(hx, px)  

  px2 = seq.convolute(hx, px1)

  py1 = seq.convolute(hy, py)  

  py2 = seq.convolute(hy, py1)

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


newton.foldchange.MK = function (x, y, Nx, Ny, w1 = 1.0, w2 = 1.0, tol = 1e-09, init.base = 1.5e-8, init.fold = 1.6, 
                                         max.iter = 30, no.log=TRUE, show.iter = FALSE) {   

  LARGE_M = 500; SMALL_M = 1e-9  # to avoid over/under flow, Nov 8, 2020

  beta0 = matrix(c(log(init.base), log(init.fold)), ncol = 2)

  if (show.iter) 

    message(paste("iteration 0 yielding ... ", toString(beta0)))

  for (i in 1:max.iter) {

  if ((Nx*exp(beta0[1]) > LARGE_M) || (Ny*exp(beta0[1] + beta0[2]) > LARGE_M) || 
      (Nx*exp(beta0[1]) < SMALL_M) || (Ny*exp(beta0[1] + beta0[2]) < SMALL_M))  {

     message('A very large/small value of m has been encountered. New starting value may help.'); return(NA)}  # Nov 8, 202

  srinfo = MK.fold.score.info(x, y, Nx, Ny, beta0[1], beta0[2], w1 = w1, w2 = w2)

  score = matrix(c(srinfo[1], srinfo[2]), ncol=2) 

  fisher = matrix(c(srinfo[3], srinfo[4], srinfo[4], srinfo[5]), ncol = 2, byrow = TRUE)   

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

}  # end of newton.foldchange.MK



### ========================== profile CI for b1, mutation rate fold change, Nov 3, 2020 ===========

# ----------- define starting value for constructing CI for b1, Nov 3, 2010 ----------

fold.MK.venzon.b1 = function(x, y, Nx, Ny, b0.hat, b1.hat, alpha=0.05, w1 = 1.0, w2 = 1.0) {

  srinfo = MK.fold.score.info(x, y, Nx, Ny, b0.hat, b1.hat, w1 = w1, w2 = w2)

  q = qchisq(1 - alpha, df = 1)

  J11 = srinfo[3]; J12 = srinfo[4]; J22 = srinfo[5]

  h = J22 - J12^2 / J11

  h = 0.5 * sqrt(q / h)

  ven.low.b0 = b0.hat + h * J12 / J11; ven.low.b1 = b1.hat - h

  ven.up.b0 = b0.hat - h * J12 / J11; ven.up.b1 = b1.hat + h

return( c(ven.low.b0, ven.low.b1, ven.up.b0, ven.up.b1) ) 
} 


# ----------- the log likelihood needed for calculating init values, Nov 3, 2020 ---------------

fold.loglik.MK = function(x, y, Nx, Ny, b0, b1, w1 = 0.8, w2 = 0.8) {

  mx = Nx*exp(b0)

  my = Ny * exp(b0 + b1)

  probx = prob.MK(m = mx, w = w1, n = max(x))

  proby = prob.MK(m = my, w = w2, n = max(y))

  loglik = sum(log(probx[x + 1])) + sum(log(proby[y + 1]))

return(loglik)

}


# --------------- computing profile CI for b1, the logarithm of fold change, Nov 3, 2020 ---------------

confint.foldchange.MK = function(x, y, Nx, Ny, w1 = 1.0, w2 = 1.0, alpha = 0.05, init.base = 1e-8,
                                 init.fold = 1.6, init.low.base = -9, init.up.base = -9, init.low.fold = -9,
                                 init.up.fold = -9, max.iter = 30, tol = 1e-9, show.iter = FALSE) {

  init.b00 = log(init.base)

  init.b10 = log(init.fold)

  if (show.iter) cat('\nComputing MLE ...\n')

  mle = newton.foldchange.MK(x, y, Nx, Ny, w1 = w1, w2 = w2, init.base = init.base, init.fold = init.fold, 
                                 tol = tol, max.iter = max.iter, no.log = FALSE, show.iter = show.iter)

  if (is.na(mle)[1]) return(NA) # Nov 8, 2020

  if (show.iter) cat('MLE of b0 and b1 ...',  mle, ' \n')

  b0.hat = mle[1]; b1.hat = mle[2]

  if ((init.low.base < 0) || (init.low.fold < 0) || (init.up.base < 0) || (init.up.base < 0))

    {all.inits = fold.MK.venzon.b1(x, y, Nx, Ny, b0.hat, b1.hat, w1 = w1, w2 = w2, alpha = alpha) }

  if (init.low.base < 0) {b00 = all.inits[1]} else b00 = log(init.low.base)

  if (init.low.fold < 0) {b10 = all.inits[2]} else b10 = log(init.low.fold)


### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

  srinfo = MK.fold.score.info(x, y, Nx, Ny, b00, b10, w1 = w1, w2 = w2)

  zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

  like = fold.loglik.MK(x, y, Nx, Ny, b00, b10, w1 = w1, w2 = w2)

  la = fold.loglik.MK(x, y, Nx, Ny, b0.hat, b1.hat, w1 = w1, w2 = w2) - 0.5 * qchisq(1 - alpha, df = 1)

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

    srinfo = MK.fold.score.info(x, y, Nx, Ny, b01, b11, w1 = w1, w2 = w2)

    zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

    like=fold.loglik.MK(x, y, Nx, Ny, b01, b11, w1 = w1, w2 = w2)

    zheng37vet=matrix(c(like - la, srinfo[1]), ncol = 1)

    beta0 = beta1

    } ## end of for lower limit

# --------------- repeating for the upper limit --------------------


  if (init.up.base<0) {b00 = all.inits[3]} else b00 = log(init.up.base)

  if (init.up.fold<0) {b10 = all.inits[4]} else b10 = log(init.up.fold)

### The matrix and vector in equation 37 of Zheng (Math. Biosci, 2015) 

  srinfo = MK.fold.score.info(x, y, Nx, Ny, b00, b10, w1 = w1, w2 = w2)

  zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

  like = fold.loglik.MK(x, y, Nx, Ny, b00, b10, w1 = w1, w2 = w2)

  la = fold.loglik.MK(x, y, Nx, Ny, b0.hat, b1.hat, w1 = w1, w2 = w2) - 0.5 * qchisq(1 - alpha, df = 1)

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

    srinfo = MK.fold.score.info(x, y, Nx, Ny, b01, b11, w1 = w1, w2 = w2)

    zheng37mat = matrix(c(-srinfo[1], -srinfo[2], srinfo[3], srinfo[4]), ncol = 2, byrow = TRUE)

    like = fold.loglik.MK(x, y, Nx, Ny, b01, b11, w1 = w1, w2 = w2)

    zheng37vet = matrix(c(like - la, srinfo[1]), ncol = 1)

    beta0 = beta1

    } ## end of for upper limit

  if (show.iter) cat('\n')

  return(exp(c(b1.lower, mle[2], b1.upper)))

}  ########## end of confint.foldchange.LD, Nov 4, 2020


