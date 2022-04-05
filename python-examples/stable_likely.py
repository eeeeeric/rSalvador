# Aug 19, 2021, starting maximum likelihoood estimation
# Feb 17, 2022, using faster algorithmic approach to computing the probability function

import numpy as np
import mpmath as mp

################ computing the probability function ###############

def find_g_seq(n=8, m=2.5, dr=0.2):
   h=(n+1)*[0]
   g=(n+1)*[0]
   h[1]=1
   g[1]=m/(2-dr)
   for k in range(2,n+1):
      h[k]=h[k-1]*(k-1)/(k-dr)
      g[k]=h[k]*m/(k+1-dr)
   return np.array(g) 


# Feb 16, 2022: faster method 
def probKL_stable (k=3, m=2.1, dr=0.3):
   if (k<0):
      return None
   else:
      if (k==0):
         return [np.exp(-m/(1-dr))]
   prob=(k+1)*[0]
   prob[0]=np.exp(-m/(1-dr))
   prob[1]=m*np.exp(-m/(1-dr))/(2-dr)
   g=find_g_seq(k, m, dr)
   for n in range(2,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*g[j]*prob[n-j]
      prob[n]=prob[n]/n
   return np.array(prob)

############### end of probability computing ##################


# Aug 21, 2021, simplify the model

def loglikely(Nt,Y,mu,dr):
   ntube=len(Nt)
   loglik=0
   for i in range(ntube):
      m=Nt[i]*mu
      prob=probKL_stable(k=Y[i],m=m,dr=dr)[-1]
      loglik=loglik+mp.log(prob)
   return(float(loglik))



def golden_KL0(data, Nt, dr=0.2, mu_low=1e-10, mu_up=1e-4, tol=1e-5, max_iter=75, show_iter=True):
   gr=(np.sqrt(5)+1)/2
   a=mu_low
   b=mu_up
   c = b - (b - a)/gr
   d = a + (b - a)/gr
   if show_iter:
      print([a,b])
   for i in range(max_iter):
      fc= loglikely(Nt=Nt,Y=data,mu=c,dr=dr)
      fd= loglikely(Nt=Nt,Y=data,mu=d,dr=dr)
      if (fd < fc):
         b=d
      else:
         a=c
      c=b-(b-a)/gr
      d=a+(b-a)/gr
      if show_iter:
         print('iteration '+str(i)+' yields ... '+str([a,b]) )
      if abs(b-a)<tol*(abs(c)+abs(d)):
         return( (a+b)/2 )
   return(-9.999) #when it exceed max number of iteration



