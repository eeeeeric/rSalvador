# Apr 12, 2022, bisection to get confidence intervals for the plating+fitness model
# Jul 14, 2022, let Nt, e and w be colony-dependent, to analyze Gasch's data
import scipy.stats as pystats

import gasch_likely as lik  

def confint_bisect(data, Nt=[123], w=[0.9], e=[0.1], alpha=0.05,\
     mu_low=1e-8, mu_up=0.1, low_fac=0.1, up_fac=1.8, tol=1e-3,\
     max_iter=90, show_iter=True):
   qa=pystats.chi2.ppf(1-alpha,1)
   mu_hat=lik.golden_plafit(data, Nt=Nt, w=w, e=e, mu_low=mu_low, mu_up=mu_up, tol=tol)
   delta=-lik.plafit_loglikely(data, Nt,  w, e, mu_hat) + qa/2

   def f(mu):
      return(lik.plafit_loglikely(data, Nt, w, e, mu)+delta)
      #  end of inner function f

   left=low_fac*mu_hat
   right=mu_hat

   for i in range(max_iter):
      c=(left+right)/2
      fc=f(c)
      if (f(c)<0): left=c
      else: right=c
      if show_iter: print('iteration '+str(i)+' yielding '+str(float(c))+ ' ...')
      if (abs(left-right)/mu_hat)<tol: low_lim=(left+right)/2; break

# now the upper CI limit
   right=up_fac*mu_hat
   left=mu_hat

   for i in range(max_iter):
      c=(left+right)/2
      fc=f(c)
      if (f(c)<0): right=c
      else: left=c
      if show_iter: print('iteration '+str(i)+' yielding '+str(float(c))+ ' ...')
      if (abs(left-right)/mu_hat)<tol: up_lim=(left+right)/2; break

   return([float(mu_hat), float(low_lim), float(up_lim)])


