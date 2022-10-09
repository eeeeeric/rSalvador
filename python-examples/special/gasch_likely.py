# Apr 17, 2022, stable likelihood for the model with fitness differential and imperfect plating efficiency
# Jul 14, 2022, let Nt, e, and w vary across colonies to accommodate Gasch's data
import numpy as np
import mpmath as mp
mp.mp.dps = 40

import simplest as Bell

#Apr 12, 2022

def plafit_loglikely(Y,Nt,w,e,mu):
   ncolony=len(Y)
   prob=np.array([-9.9]*ncolony)
   for i in range(ncolony):
      prob[i]=Bell.prob_plafit(k=Y[i], m=Nt[i]*mu, w=w[i], e=e[i])[-1]
      prob[i]=np.log(prob[i])
   return(np.sum(prob))


def golden_plafit(data, Nt=[1234], w=[0.75], e=[0.1], mu_low=1e-9, mu_up=0.1, tol=1e-5, max_iter=75, show_iter=True):
   gr=(mp.sqrt(5)+1)/2
   a=mu_low
   b=mu_up
   c = b - (b - a)/gr
   d = a + (b - a)/gr
   if show_iter:
      print([a,b])
   for i in range(max_iter):
      fc= plafit_loglikely(Y=data, Nt=Nt, w=w, e=e, mu=c)
      fd= plafit_loglikely(Y=data, Nt=Nt, w=w, e=e, mu=d)
      if (fd < fc):
         b=d
      else:
         a=c
      c=b-(b-a)/gr
      d=a+(b-a)/gr
      if show_iter:
         print('iteration '+str(i)+' yields ... '+str([float(a),float(b)]) )
      if abs(b-a)<tol*(abs(c)+abs(d)):
         return( (a+b)/2 )
   return(-9.999) #when it exceed max number of iteration

