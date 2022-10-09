# starting Aug 11, 2021

import mpmath as mp

#mp.mp.dps = 20
mp.mp.dps = 100
#mp.mp.dps = 5275


# integrand for the case where both fitness and plating efficiency must be considered

def plat_integrand(t=1.1, r=0.8, k=2, m=2.2, w=0.9, e=0.1):
   exp_it=complex(mp.cos(t), mp.sin(t))
   exp_ikt=complex(mp.cos(k*t), mp.sin(k*t))
   A=mp.hyp2f1(1, 1/w, 1+1/w, ( (1-e)/e +r*exp_it)/(-1+r*exp_it))
   return mp.exp(-m*A)/r**k/exp_ikt


# computing LD probability by integrating the KL integrand

def prob_cauchy(k=2, m=2.2, w=0.9, e=0.1, r=0.9):

   prob= mp.quad(lambda t:plat_integrand(t=t, r=r, k=k, m=m, w=w, e=e), [0, 2*mp.pi])
   return(float(prob.real/2/mp.pi))

