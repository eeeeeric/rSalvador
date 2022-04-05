# starting Aug 11, 2021

import mpmath as mp


mp.mp.dps = 40
#mp.mp.dps = 50
#mp.mp.dps = 5275

# Nov 15, 2021, for the example probLD_kl(0,0.008), is hard, and hence need large dps

# integrand given by Kessler and Levine

def KL_luria_integrand(t,m,k):

   return (1+t)**(-k-1) *t**(-m*t/(1+t)) *mp.sin(mp.pi*m*t/(1+t))/mp.pi


# computing LD probability by integrating the KL integrand

def probLD_kl(k=6,m=2.3):

   prob= mp.quad(lambda t:KL_luria_integrand(t,m,k), [0, mp.inf])

   return(float(prob))


def KL_koch_integrand(t,k,m,r):
   a1=m*r*t/(r-1)/(1+t) *mp.hyp2f1(1,1-r,2-r,t/(1+t))
   a2=m*mp.pi*r/mp.tan(mp.pi*r)*(t/(1+t))**r
   a3=m*mp.pi*r*(t/(1+t))**r
   a=mp.exp(a1)*mp.exp(-a2)*mp.sin(a3)/(1+t)**(k+1) / mp.pi
   return a


def probMK_kl(k=6,m=1.3,r=1.2):

   ''' probMK_kl(k=6,m=1.3,r=1.2) computes mutant probabilities using equation (29)
   of Kessler and Levine (2015). Here r is the ratio of birth rate of wild-type
   cells to that of mutant cells, which is the reciprocal of fitness. '''

   prob= mp.quad(lambda t:KL_koch_integrand(t, k, m, r), [0, mp.inf])
   return(float(prob))


# Page 792 of Kessler and Levine (2015), eq (48), bw is birth rate of wile-type, etc

def KL_koch_death_integrand(t,k,m,bw,dw,bm,dm):
   r=(bw-dw)/(bm-dm)
   a1=m*bw*t/(r-1)/(bm*(1+t)-dm) *mp.hyp2f1(1,1-r,2-r,bm*t/(bm*(1+t)-dm))
   a2=m*mp.pi*bw/bm/mp.tan(mp.pi*r)*(bm*t/(bm*(1+t)-dm))**r
   a3=m*mp.pi*bw/bm*(bm*t/(bm*(1+t)-dm))**r
   a=mp.exp(a1)*mp.exp(-a2)*mp.sin(a3)/(1+t)**(k+1) / mp.pi
   return a

   
def probMK_death(k=6,m=1.3,bw=1.2,dw=0.2,bm=1.2,dm=0.3):

   ''' probMK_death(k=6,m=1.3,bw=1.2,dw=0.2,bm=1.2,dm=0.3) computes mutant probabilities
   using equation (48) of Kessler and Levine (2015). bw is birth rate of wild-type cells,
   dw is death rate of wild-type cells, bw is birth rate of mutants, and dm is death rate
   of mutants. '''

   prob= mp.quad(lambda t:KL_koch_death_integrand(t, k, m, bw,dw,bm,dm), [0, mp.inf])
   return(float(prob))


# Aug 21, 2021: simplified KL nodel for subinhibitory exposure

def simple_KL_integrand(t,k,m,dr):
   a1=m*t/dr/(1+t) *mp.hyp2f1(1,dr,1+dr,t/(1+t))
   a2=m*mp.pi/mp.tan(mp.pi*(1-dr))*(t/(1+t))**(1-dr)
   a3=m*mp.pi*(t/(1+t))**(1-dr)
   a=mp.exp(-a1-a2)*mp.sin(a3)/(1+t)**(k+1) / mp.pi  ### simplified Feb 15, 2022
   return a

   
def probKL_simple(k=6,m=1.3,dr=0.2):

   ''' probKL_simple(k=6,m=1.3,dr=0.2) computes mutant probabilities
   using equation (48) of Kessler and Levine (2015). dr is the relative death rate
   defined by dw/bw. '''

   prob= mp.quad(lambda t:simple_KL_integrand(t, k, m, dr), [0, mp.inf])
   return(float(prob))




