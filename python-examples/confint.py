# Aug 30, 2021, bisection to get confidence intervals.

import stable_likely as lik

def confint_bisect(data, Nt, dr=0.2, mu_low=1e-8, mu_up=9e-8,\
    low_fac=0.1, up_fac=1.8, tol=1e-3, max_iter=90, show_iter=True):

   mu_hat=lik.golden_KL0(data, Nt, dr=dr, mu_low=mu_low, mu_up=mu_up, tol=tol)
   delta=-lik.loglikely(Nt, data, mu_hat, dr)+1.9207

   def f(mu):
      return(lik.loglikely(Nt, data, mu, dr)+delta)
      #  end of inner function f

   left=low_fac*mu_hat
   right=mu_hat

   for i in range(max_iter):
      c=(left+right)/2
      fc=f(c)
      if (f(c)<0): left=c
      else: right=c
      if show_iter: print('iteration '+str(i)+' yielding '+str(c)+ ' ...')
      if (abs(left-right)/mu_hat)<tol: low_lim=(left+right)/2; break

# now the upper CI limit
   right=up_fac*mu_hat
   left=mu_hat

   for i in range(max_iter):
      c=(left+right)/2
      fc=f(c)
      if (f(c)<0): right=c
      else: left=c
      if show_iter: print('iteration '+str(i)+' yielding '+str(c)+ ' ...')
      if (abs(left-right)/mu_hat)<tol: up_lim=(left+right)/2; break

   return([mu_hat, low_lim, up_lim])









 




