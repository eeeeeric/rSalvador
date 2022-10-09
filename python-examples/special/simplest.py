# Apr 17, 2022, improving results from yesterday, the most practical approach
# Apr 18, 2022: simplest code
# Jul 6, 2022: further simplification

import mpmath as mp
mp.mp.dps = 40

def eta_seq(e=0.1, w=0.9, n=3):
   eta=[0]*(n+1)
   eta[1]=mp.mpf(1)
   for k in range(1,n):
     eta[k+1]=k*e/(k+1/w) *eta[k] 
   return eta

# the f sequence
def f_seq(e=0.1, w=0.9, n=6):
   f=[mp.hyp2f1(k+1, k+1, k+1+1/w, 1-e) for k in range(1,n+1)] 
   return f

# the q sequence
def q_seq(m=2.1, e=0.1, w=0.9, n=6):
   f0=mp.hyp2f1(1,1,1+1/w, 1-e)
   f=f_seq(e=e, w=w, n=n+1)
   f.insert(0,f0)
   eta=eta_seq(e=e, w=w, n=n+1)
   q=[0]*(n+1)
   q[0]=-m*e*f0
   for k in range(1,n+1):
      q[k]=m*e*eta[k]*(f[k-1]-f[k]*k*e/(k+1/w))
   return q


# Apr 16, 2022: old-school approach

def prob_plafit (k=3, m=2.1, w=0.9, e=0.1):
   if (k<0):
      return None
   else:
      if (k==0):
         return [mp.exp(q_seq(m=m, w=w, e=e,n=0)[-1])]
   theta=(1-e)/e
   q=q_seq(m=m, e=e, w=w, n=k+1)
   prob=(k+1)*[0]
   prob[0]=mp.exp(q[0])
   for n in range(1,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*q[j]*prob[n-j]
      prob[n]=prob[n]/n
   prob=map(lambda u:float(u), prob)
   return list(prob)

