## Apr 16, 2000, using Bruno formula to handle fitness+plating

import numpy as np
import scipy.special as sc
import sympy.functions.combinatorial.numbers as combo
import math

# the f sequence
def f_seq(e=0.1, n=6):
   theta=(1-e)/e
   f=[-(1+theta)*math.factorial(k) for k in range(1,n+1)]
   return f

# the g sequence
def g_seq(e=0.1, w=0.9, n=6):
   theta=(1-e)/e
   g=[sc.hyp2f1(k+1, k+1/w, k+1+1/w, -theta)/( k*w +1) *math.factorial(k) for k in range(1,n+1)] 
   return g


# the h sequence
def h_seq(e=0.1, w=0.9, n=6):
   f=f_seq(e=e, n=n)
   g=g_seq(e=e, w=w, n=n)
   h=[0]*n
   for k in range(0,n):
      for r in range(0,k+1):
         h[k]=h[k]+combo.bell(k+1, r+1, f[0:k+1])*g[r]
   return h


# Apr 16, 2022: old-school approach

def probBruno(k=3, m=2.1, w=0.9, e=0.1):
   if (k<0):
      return None
   else:
      if (k==0):
         return [123456]
   theta=(1-e)/e
   q0=-m*sc.hyp2f1(1,1/w,1+1/w,-theta)
   q=h_seq(e=e, w=w, n=k)
   q.insert(0, q0)
   div=[(-1/m)*math.factorial(r) for r in range(0,k+1)]
   q=np.divide(q, div) 
   prob=(k+1)*[0]
   prob[0]=math.exp(q0)
   for n in range(1,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*q[j]*prob[n-j]
      prob[n]=prob[n]/n
   return np.array(prob)

