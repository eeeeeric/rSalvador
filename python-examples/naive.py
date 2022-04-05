# starting Aug 11, 2021
# Feb 13, 2022, deleting arbitrary-precision computation for comparison purposes

import numpy as np
import scipy.integrate as integrate
import scipy.special as sc

# Nov 15, 2021
# integrand given by Kessler and Levine

def KL_luria_integrand(t,m,k):
   return (1+t)**(-k-1) *t**(-m*t/(1+t)) *np.sin(np.pi*m*t/(1+t))/np.pi

# computing LD probability by integrating the KL integrand

def probLD_kl(k=6,m=2.3):
   prob= integrate.quad(lambda t:KL_luria_integrand(t,m,k), 0, np.inf)
   return(prob[0])

def KL_koch_integrand(t,k,m,r):
   a1=m*r*t/(r-1)/(1+t) *sc.hyp2f1(1,1-r,2-r,t/(1+t))
   a2=m*np.pi*r/np.tan(np.pi*r)*(t/(1+t))**r
   a3=m*np.pi*r*(t/(1+t))**r
   a=np.exp(a1)*np.exp(-a2)*np.sin(a3)/(1+t)**(k+1) / np.pi
   return a

def probMK_kl(k=6,m=1.3,r=1.2):

   ''' probMK_kl(k=6,m=1.3,r=1.2) computes mutant probabilities using equation (29)
   of Kessler and Levine (2015). Here r is the ratio of birth rate of wild-type
   cells to that of mutant cells, which is the reciprocal of fitness. '''

   prob= integrate.quad(lambda t:KL_koch_integrand(t, k, m, r), 0, np.inf)
   return(prob[0])


def probMK_kl(k=6,m=1.3,r=1.2):

   ''' probMK_kl(k=6,m=1.3,r=1.2) computes mutant probabilities using equation (29)
   of Kessler and Levine (2015). Here r is the ratio of birth rate of wild-type
   cells to that of mutant cells, which is the reciprocal of fitness. '''

   prob= integrate.quad(lambda t:KL_koch_integrand(t, k, m, r), 0, np.inf)
   return(prob[0])


# Page 792 of Kessler and Levine (2015), eq (48), bw is birth rate of wile-type, etc

def KL_koch_death_integrand(t,k,m,bw,dw,bm,dm):
   r=(bw-dw)/(bm-dm)
   a1=m*bw*t/(r-1)/(bm*(1+t)-dm) *sc.hyp2f1(1,1-r,2-r,bm*t/(bm*(1+t)-dm))
   a2=m*np.pi*bw/bm/np.tan(np.pi*r)*(bm*t/(bm*(1+t)-dm))**r
   a3=m*np.pi*bw/bm*(bm*t/(bm*(1+t)-dm))**r
   a=np.exp(a1)*np.exp(-a2)*np.sin(a3)/(1+t)**(k+1) / np.pi
   return a

   
def probMK_death(k=6,m=1.3,bw=1.2,dw=0.2,bm=1.2,dm=0.3):

   ''' probMK_death(k=6,m=1.3,bw=1.2,dw=0.2,bm=1.2,dm=0.3) computes mutant probabilities
   using equation (48) of Kessler and Levine (2015). bw is birth rate of wild-type cells,
   dw is death rate of wild-type cells, bw is birth rate of mutants, and dm is death rate
   of mutants. '''

   prob= integrate.quad(lambda t:KL_koch_death_integrand(t, k, m, bw,dw,bm,dm), 0, np.inf)
   return(prob[0])


# Aug 21, 2021: simplified KL nodel for subinhibitory exposure

def simple_KL_integrand(t,k,m,dr):
   a1=m*t/dr/(1+t) *sc.hyp2f1(1,dr,1+dr,t/(1+t))
   a2=m*np.pi/np.tan(np.pi*(1-dr))*(t/(1+t))**(1-dr)
   a3=m*np.pi*(t/(1+t))**(1-dr)
   a=np.exp(-a1-a2)*np.sin(a3)/(1+t)**(k+1) / np.pi  ## simplified Feb 15, 2022
   return a

   
def probKL_simple(k=6,m=1.3,dr=0.2):

   ''' probKL_simple(k=6,m=1.3,dr=0.2) computes mutant probabilities
   using equation (48) of Kessler and Levine (2015). dr is the relative death rate
   defined by dw/bw. '''

   prob= integrate.quad(lambda t:simple_KL_integrand(t, k, m, dr), 0, np.inf)
   return(prob[0])





