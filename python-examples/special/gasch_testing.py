# Jul 17, 2022, likelihood ratio test for Gasch's data

import numpy as np

import scipy.stats as pystats

import gasch_likely as Gold

def plafit_LR(Y1, Y2, Nt1, Nt2, w1, w2, e1, e2,\
              mu_low1=1e-9, mu_low2=1e-9, mu_lowc=1e-9,\
              mu_up1=0.002, mu_up2=0.002, mu_upc=0.002):
   Y=list(Y1)+list(Y2)
   Nt=list(Nt1)+list(Nt2)
   w=list(w1)+list(w2)
   e=list(e1)+list(e2)

   hat_mu1=Gold.golden_plafit(Y1, Nt1, w1, e1, mu_low=mu_low1, mu_up=mu_up1)
   hat_mu2=Gold.golden_plafit(Y2, Nt2, w2, e2, mu_low=mu_low2, mu_up=mu_up2)
   hat_muc=Gold.golden_plafit(Y, Nt, w, e, mu_low=mu_lowc, mu_up=mu_upc)

   lik1=Gold.plafit_loglikely(Y1, Nt1, w1, e1, hat_mu1)
   lik2=Gold.plafit_loglikely(Y2, Nt2, w2, e2, hat_mu2)
   likc=Gold.plafit_loglikely(Y, Nt, w, e, hat_muc)
   lam=2*(lik1+lik2-likc)
   pval=pystats.chi2.sf(x=lam, df=1)
   return([lam, pval])





   
