# Feb 15, 2022, first few terms of the PGF

import scipy.special as sp
import numpy as np
import math

def KL_p0(m=2.5, dr=0.2):
   return np.exp(-m/(1-dr))

def KL_p1(m=2.5, dr=0.2):
   return np.exp(-m/(1-dr))*m/(2-dr)

def KL_p2(m=2.6, dr=0.2):
   return 0.5*np.exp(-m/(1-dr))*(m**2/(2-dr)**2 +2*m/(2-dr)/(3-dr) )

def KL_p3(m=2.6, dr=0.2):
   a=48+48*m+12*m**2-48*dr-36*m*dr-7*m**2*dr+12*dr**2+6*m*dr**2+m**2*dr**2
   b=(dr-4)*(dr-3)*(2-dr)**3
   return np.exp(-m/(1-dr))*m*a/b/6


# Mar 6, 2022, exact probability p4, code adapted from Mathematica output
def KL_p4(m=2.6, dr=0.2):
   A=3456 - 6336*dr + 4320*dr**2 - 1296*dr**3 + \
      144*dr**4 + 3840*m - 5808*dr*m + 3168*dr**2*m - 732*dr**3*m + \
      60*dr**4*m + 1440*m**2 - 1848*dr*m**2 + 852*dr**2*m**2 - \
      168*dr**3*m**2 + 12*dr**4*m**2 + 180*m**3 - 201*dr*m**3 + \
      83*dr**2*m**3 - 15*dr**3*m**3 + dr**4*m**3
   B=24*(-5 + dr)*(-4 + dr)*(-3 + dr)**2*(-2 + dr)**4
   prob=np.exp(-m/(1-dr))*m*A/B
   return prob

# Mar 5, 2022, for p5
def KL_p5(m=2.6, dr=0.2):
   A= -138240*m + 322560*dr*m - 299520*dr**2*m + \
      138240*dr**3*m - 31680*dr**4*m + 2880*dr**5*m - 161280*m**2 + \
      314880*dr*m**2 - 238080*dr**2*m**2 + 86400*dr**3*m**2 - \
      14880*dr**4*m**2 + 960*dr**5*m**2 - 72000*m**3 + 120000*dr*m**3 - \
      76320*dr**2*m**3 + 23040*dr**3*m**3 - 3300*dr**4*m**3 + \
      180*dr**5*m**3 - 14400*m**4 + 20880*dr*m**4 - 11600*dr**2*m**4 + \
      3100*dr**3*m**4 - 400*dr**4*m**4 + 20*dr**5*m**4 - 1080*m**5 + \
      1386*dr*m**5 - 699*dr**2*m**5 + 173*dr**3*m**5 - 21*dr**4*m**5 + dr**5*m**5
   B= 120*(-6 + dr)*(-5 + dr)*(-4 + dr)*(-3 + dr)**2* (-2 + dr)**5
   prob=-np.exp(-m/(1-dr))*A/B
   return prob


def asympKL(k=123, m=2.5, dr=0.2):
   return m*math.gamma(2-dr)/k**(2-dr)


# temporary, Feb 16, 2022
def coef_g(k=3, m=2.5, dr=0.2):
   return math.factorial(k-1)*m/(k+1-dr)/sp.poch(2-dr,k-1)

def g_seq(n=8, m=2.5, dr=0.2):
   B=(n+1)*[0]
   B[1]=m/(2-dr)
   for k in range(2,n+1):
      B[k]=coef_g(k, m, dr)
   return np.array(B) 


def find_g_seq(n=8, m=2.5, dr=0.2):
   h=(n+1)*[0]
   g=(n+1)*[0]
   h[1]=1
   g[1]=m/(2-dr)
   for k in range(2,n+1):
      h[k]=h[k-1]*(k-1)/(k-dr)
      g[k]=h[k]*m/(k+1-dr)
   return np.array(g) 



# Feb 16, 2022: probability functin of the Kessler distribution
def probKL (k=3, m=2.1, dr=0.3):
   if (k<0):
      return None
   else:
      if (k==0):
         return np.exp(-m/(1-dr))
   prob=(k+1)*[0]
   prob[0]=np.exp(-m/(1-dr))
   prob[1]=m*np.exp(-m/(1-dr))/(2-dr)
   B=g_seq(k, m, dr)
   for n in range(2,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*B[j]*prob[n-j]
      prob[n]=prob[n]/n
   return np.array(prob)


# Feb 16, 2022: faster method 
def slapdash (k=3, m=2.1, dr=0.3):
   if (k<0):
      return None
   else:
      if (k==0):
         return [np.exp(-m/(1-dr))]
   prob=(k+1)*[0]
   prob[0]=np.exp(-m/(1-dr))
   prob[1]=m*np.exp(-m/(1-dr))/(2-dr)
   B=find_g_seq(k, m, dr)
   for n in range(2,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*B[j]*prob[n-j]
      prob[n]=prob[n]/n
   return np.array(prob)

