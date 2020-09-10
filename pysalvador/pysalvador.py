# pySalvador 0.1, February 14, 2020, Qi Zheng, Texas A&M University School of Public Health, College Station, TX 77843

import numpy as np
import scipy.special as sc
import scipy.stats as pystats


# Feb 6, 2020: CI for the MK model
def confintMK (data,w=1,alpha=0.05,tol=1e-8,init_m=0,init_lower=0,init_upper=0,max_iter=30,show_iter=False):
   initLow=init_lower
   initUp=init_upper
   hseq=betaSeq(w,max(data))/w
   hseq[0]=-1
   if (show_iter):
      print('Iterating for MLE of m ... ')
   mhat=newtonMK(data,w=w,tol=tol,init_m=init_m,show_iter=show_iter)
   if (show_iter):
      print('ML estimate of m is ... '+str(mhat))
   score_info_like=MK_score_info_loglikely(mhat,w=w,hseq=hseq,data=data)
   score=score_info_like[0]
   info=score_info_like[1]
   like=score_info_like[2]
   qa=pystats.chi2.ppf(1-alpha,1)
   h=np.sqrt(qa/info)
   la=like-0.5*qa
   if (show_iter):
      print('Iterating for lower limit ... ')
   if (initLow<=0):
      m0=mhat-0.5*h
   else:
      m0=initLow
   for i in range(1,max_iter+1):
      like_score=MK_loglike_score(m0,w,hseq,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs((m1-m0)/m0)<tol):
         mL=m1
         break
      else:
         m0=m1
# now upper limit
   if (show_iter):
      print('Iterating for upper limit ... ')
   if (initUp<=0):
      m0=mhat+0.5*h
   else:
      m0=initUp
   for i in range(1,max_iter+1):
      like_score=MK_loglike_score(m0,w,hseq,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs(m1-m0)/m0<tol):
         mU=m1
         break
      else:
         m0=m1
   return [mL,mU]
      

# Feb 6, 2020: score, info and likelihood for Mandelbrot-Koch model
def MK_score_info_loglikely (m,w,hseq,data): # keep hseq out of looping
   n=max(data)
   p=probMK(m,w,n)
   p1=seq_convolute(hseq,p)
   p2=seq_convolute(hseq,p1)
   loglike=np.log([p[i] for i in data])
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info),sum(loglike)]


# Feb 6, 2020: reduced version of MK_score_info_loglikely
def MK_score_info (m,w,hseq,data): # keep hseq out of looping
   n=max(data)
   p=probMK(m,w,n)
   p1=seq_convolute(hseq,p)
   p2=seq_convolute(hseq,p1)
   score=[(p1/p)[i] for i in data]
   info=[((p1/1)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info)]


def MK_loglike_score (m,w,hseq,data): # keep hseq out of looping
   n=max(data)
   p=probMK(m,w,n)
   loglike=np.log([p[i] for i in data])
   p1=seq_convolute(hseq,p)
   score=[(p1/p)[i] for i in data]
   return [sum(loglike),sum(score)]


# Feb 6, 2020: likelihodd and score for the Mandelbrot-Koch model
def MK_loglike_score999 (m,w,hseq,data): # keep hseq out of looping
   n=max(data)
   p=probMK(m,w,n)
   p1=seq_convolute(hseq,p)
   loglike=np.log([p[i] for i in data])
   score=[(p1/p)[i] for i in data]
   return [sum(loglike),sum(score)]


# Feb 5, 2020: CI for LD with parting plating
def confintLD_plating (data,e,alpha=0.05,tol=1e-6,init_m=0,init_lower=0,init_upper=0,max_iter=30,show_iter=False):
   initLow=init_lower
   initUp=init_upper
   eta=etaSeq(e,max(data))
   if (show_iter):
      print('Iterating for MLE of m ... ')
   mhat=newtonLD_plating(data,e=e,tol=tol,init_m=init_m,show_iter=show_iter)
   if (show_iter):
      print('The ML estimate of m ... '+str(mhat))
   score_info_like=LD_score_info_loglikely_plating(mhat,e=e,eta=eta,data=data)
   score=score_info_like[0]
   info=score_info_like[1]
   like=score_info_like[2]
   qa=pystats.chi2.ppf(1-alpha,1)
   h=np.sqrt(qa/info)
   la=like-0.5*qa
   if (show_iter):
      print('Iterating for lower limit ... ')
   if (initLow<=0):
      m0=mhat-0.5*h
   else:
      m0=initLow
   for i in range(1,max_iter+1):
      like_score=LD_loglike_score_plating(m0,e,eta,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs(m1-m0)/m0<tol):
         mL=m1
         break
      else:
         m0=m1
# end of lower limit
   if (show_iter):
      print('Iterating for upper limit ... ')
   if (initUp<=0):
      m0=mhat+0.5*h
   else:
      m0=initUp
   for i in range(1,max_iter+1):
      like_score=LD_loglike_score_plating(m0,e,eta,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs(m1-m0)/m0<tol):
         mU=m1
         break
      else:
         m0=m1
   return [mL,mU]


# Feb 5, 2020: eta needs computing only once, hence is outside this function
def LD_score_info_loglikely_plating (m,e,eta,data):
   k=e/(1-e)
   n=max(data)
   p=probLD_plating(m,e,n)
   p1=seq_convolute(eta,p)*k
   p2=seq_convolute(eta,p1)*k
   loglike=np.log([p[i] for i in data])
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info),sum(loglike)]


# Feb 5, 2020: a reduced version of LD_score_info_loglikely_plating
def LD_loglike_score_plating (m,e,eta,data):
   k=e/(1-e)
   n=max(data)
   p=probLD_plating(m,e,n)
   p1=seq_convolute(eta,p)*k
   loglike=np.log([p[i] for i in data])
   score=[(p1/p)[i] for i in data]
   return [sum(loglike),sum(score)]


# Feb 4, 2020: CI for LD model
def confintLD (data,alpha=0.05,phi=1,tol=1e-6,init_m=0,init_lower=0,init_upper=0,max_iter=30,show_iter=False):
   initLow=init_lower
   initUp=init_upper
   if (show_iter):
      print('Iterating for MLE of m ... ')
   mhat=newtonLD(data,phi=phi,tol=tol,init_m=init_m,show_iter=show_iter)
   if (show_iter):
      print('The ML estimate of M ... '+str(mhat))
   score_info_like=LD_score_info_loglikely(mhat,phi=phi,data=data)
   score=score_info_like[0]
   info=score_info_like[1]
   like=score_info_like[2]
   qa=pystats.chi2.ppf(1-alpha,1)
   h=np.sqrt(qa/info)
   la=like-0.5*qa
   if (show_iter):
      print('Iterating for lower limit ... ')
   if (initLow<=0):
      m0=mhat-0.5*h
   else:
      m0=initLow
   for i in range(1,max_iter+1):
      like_score=LD_loglike_score(m0,phi,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs((m1-m0)/m0)<tol):
         mL=m1
         break
      else:
         m0=m1
# do the same for the lower limit
   if (show_iter):
      print('Iterating for upper limit ... ')
   if (initUp<=0):
      m0=mhat+0.5*h
   else:
      m0=initUp
   for i in range(1,max_iter+1):
      like_score=LD_loglike_score(m0,phi,data)
      like=like_score[0]
      score=like_score[1]
      m1=m0-(like-la)/score
      if (show_iter):
         print('iteration '+str(i)+' yielding '+str(m1))
      if (abs((m1-m0)/m0)<tol):
         mU=m1
         break
      else:
         m0=m1
   return [mL,mU]



# Feb 4, 2020: putting likelyhood, score and info together for LD confidence interval
def LD_score_info_loglikely (m,phi,data):
   n=max(data)
   p=probLD(m,phi,n)
   loglike=np.log([p[i] for i in data])
   h=list(map(lambda x:phi**(x-1)*(1/x-phi/(x+1)), range(1,n+1)))
   h.insert(0,-1)
   p1=seq_convolute(p,h)
   p2=seq_convolute(p1,h)
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info),sum(loglike)]

# Feb 4, 2020: putting likelyhood and score together for LD, a sub-function of LD_score_info_loglikely
def LD_loglike_score (m,phi,data):
   n=max(data)
   p=probLD(m,phi,n)
   loglike=np.log([p[i] for i in data])
   h=list(map(lambda x:phi**(x-1)*(1/x-phi/(x+1)), range(1,n+1)))
   h.insert(0,-1)
   p1=seq_convolute(p,h)
   score=[(p1/p)[i] for i in data]
   return [sum(loglike),sum(score)]


# Feb 4, 2020: newton for MK model with fixed fitness w
def newtonMK (data,w=1,tol=1e-8,init_m=0,max_iter=30,show_iter=False):
   if (init_m>0):
      m0=init_m
   else:
      if (np.median(data)>0):
         m0=jones_median_plating(data,e=w)
      else:
         m0=LD_p0_est(data)
   if (show_iter):
      print('iteration 0 yielding ... '+str(m0))
   h=betaSeq(w,max(data))/w
   h[0]=-1
   for i in range(1,max_iter):
      score_info=MK_score_info(m0,w,h,data)
      m1=m0+score_info[0]/score_info[1]
      if (abs(m1-m0)/m0<tol):
         return m1
      else:
         if (show_iter):
            print('iteration '+str(i)+' yielding ... '+str(m1))
         m0=m1
   return 'No convergence'


# Feb 4, 2020: score and Fisher information for MK model with known fitness w
def MK_score_info (m,w,h,data):
   n=max(data)
   p=probMK(m,w,n)
   p1=seq_convolute(h,p)
   p2=seq_convolute(h,p1)
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info)]

# Feb 4, 2020: probability functin of the MK distribution
def probMK (m=2.1,w=1.3,k=10):
   if (k<0):
      return None
   else:
      if (k==0):
         return np.exp(-m)
   prob=(k+1)*[0]
   prob[0]=np.exp(-m)
   prob[1]=m*np.exp(-m)/(1+w)
   B=betaSeq(w,k)
   for n in range(2,k+1):
      for j in range(1,n+1):
         prob[n]=prob[n]+j*B[j]*prob[n-j]
      prob[n]=prob[n]*m/n/w
   return np.array(prob)


# Feb 4, 2020: the B (beta) sequence as defined on p. 205 of Math Biosci 196 (2005) 198-204
# Leaving B[0] empty for convenience
def betaSeq (w,n):
   r=1/w
   if (n==1):
      return 1/(1+r)
   B=(n+1)*[0]
   B[1]=1/(1+r)
   for k in range(2,n+1):
      B[k]=(k-1)/(k+r) * B[k-1]
   return np.array(B)


# Feb 2, 2020: newton algorithm for LD with partial plating
def newtonLD_plating (data, e, tol=1e-08, init_m=0, max_iter=20, show_iter=False):
   if (init_m>0):
      m0=init_m
   else:
      if (np.median(data)==0):
         m0=p0_plating(data,e)
      else:
         m0=jones_median_plating(data,e)
   if (show_iter):
      print('iterating 0 yielding ... '+str(m0))
   eta=etaSeq(e,max(data))
   for i in range(1,max_iter):
      score_info=LD_score_info_plating(m0,e,eta,data)
      m1=m0+score_info[0]/score_info[1]
      if (abs(m1-m0)/m0<tol):
         return m1
      else:
         if (show_iter):
            print('iteration '+str(i)+' yielding ... '+str(m1))
         m0=m1
   return 'No convergence'


# Feb 1, 2020, P0 method for LD plating
def p0_plating (data,e):
   z=data.count(0)
   if (z==0):
      print('P0 method is not applicable.')
      return None
   else:
      n=len(data)
      est=(1-e)*(np.log(z/n)/e/np.log(e))
      return est


# Feb 1, 2020: score and info for LD with partial plating
def LD_score_info_plating (m,e,eta,data):
   k=e/(1-e)
   n=max(data)
   p=probLD_plating(m,e,n)
   p1=seq_convolute(eta,p)*k
   p2=seq_convolute(eta,p1)*k
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info)]


# Feb 1, 2020, Jones median estimator for partial plating
def jones_median_plating (data,e):
   r=np.median(data)
   if (r==0):
      print('Jones median method is not applicable.')
      return None
   est=(r/e-np.log(2))/(np.log(r/e)-np.log(np.log(2)))
   return est

# Jan 29, 2020, probability function of LD with partial plating
def probLD_plating (m=2.1, e=0.1, n=10):
   eta=etaSeq(e,n)
   odds=e/(1-e)
   prob=(n+1)*[0]
   prob[0]=np.exp(m*odds*np.log(e))
   for i in range(1,n+1):
      for j in range(1,i+1):
         prob[i]=prob[i]+j*eta[j]*prob[i-j]
      prob[i]=prob[i]*m*odds/i
   return np.array(prob)


# Jan 29, 2020, the eta sequence for partial plating
def etaSeq (e,n):
   if (e<=0.5):
      eta=etaSmall(e,n)
   else:
      eta=etaBig(e,n)
   return eta

# Jan 29, 2020, computing the eta sequence for small e<0.5, Math Biosci 216 (2008) 150-153
def etaSmall (e,n):
   odds=e/(1-e)
   eta=(n+1)*[0]
   eta[0]=np.log(e)
   eta[1]=-1-np.log(e)/(1-e)
   for i in range(1,n):
      eta[i+1]=-odds*eta[i]+1/i/(i+1)
   return np.array(eta)

# Jan 29, 2020, computing the eta sequence for small e>0.5, Math Biosci 216 (2008) 150-153
def etaBig (e,n):
   revOR=(1-e)/e
   eta=(n+1)*[0]
   eta[0]=np.log(e)
   eta[n]=((1-e)/n/(n+1))*sc.hyp2f1(1,2,n+2,1-e)
   for i in reversed(range(1,n)):
      eta[i]=revOR*(1/i/(i+1)-eta[i+1])
   return np.array(eta)

# Jan 29, 2020, simulating mutants
def simu_mutants (mu,b1,b2,N0,Nt):
   T=np.log(Nt/N0)/b1
   mean_mutation=mu*N0*(np.exp(b1*T)-1)/b1
   n_mutation=np.random.poisson(mean_mutation,1)
   mutation_epoch=np.log(1+np.random.uniform(0,1,n_mutation)*(np.exp(b1*T)-1))/b1
   offspring=map(lambda x:np.random.geometric(np.exp(-b2*(T-x))),mutation_epoch)
# numpy's definition of the hypergeometric distribution is not the same as that of R, Jan 31, 2020
   return sum(offspring)

# Jan 29, 2020, simulating sultures
def simu_cultures (n,mu,b1,b2,N0,Nt):
   cultures=n*[-1]
   for i in range(0,n):
      cultures[i]=simu_mutants(mu,b1,b2,N0,Nt)
   return cultures


# Jan 28, 2020: estimating m using an LD model
def newtonLD (data, phi=1, tol=1e-8, init_m=0, max_iter=30, show_iter=False):
   if (init_m>0):
      m0=init_m
   else:
      if (np.median(data)==0):
         m0=LD_p0_est(data)
      else:
         m0=jones_median_est(data)
   if (show_iter):
      print("iteration 0 yielding ... "+str(m0))
   for i in range(1,max_iter):
      score_info=LD_score_info(m0,phi,data)
      m1=m0+score_info[0]/score_info[1]
      if (abs(m1-m0)/m0<tol):
         return m1
      else:
         if (show_iter):
            print("iteration "+str(i)+" yielding ... "+str(m1))
         m0=m1
   return "no convergence"



# 1-28-2020, the Demerec data
demerec_data=[33,18,839,47,13,126,48,80,9,71,196,66,28,17,27,37,126,33,12,44,28,67,730,168,44,50,583,23,17,24]

# 2-13-2020, experiment 16 of Luria and Delbruck
luria_16_data=[1,0,3,0,0,5,0,5,0,6,107,0,0,0,1,0,0,64,0,35]


# Jan 28, 2020: computing probabilities of the LD distribution
def probLD (m=1.2,phi=1.0,k=9):
   problist=np.zeros(k+1)
   problist[0]=np.exp(-m)
   for n in range(1,k+1):
      cumuphi=1.0
      for j in range(1,n+1):
         problist[n]=problist[n]+cumuphi*(1-j*phi/(j+1))*problist[n-j]
         cumuphi=cumuphi*phi
      problist[n]=problist[n]*m/n
   return problist   

# Jan 28, 2020, sequence convolution
def seq_convolute (x,y):
   xlen=len(x)
   z=np.zeros(xlen)
   for k in range(1,xlen+1):
      for i in range(1,k+1):
         z[k-1]=z[k-1]+x[i-1]*y[k-i]
   return np.array(z)


# Jan 28, 2020: computing score and info matrix for the LD distribution
def LD_score_info (m, phi, data):
   n=max(data)
   p=probLD(m,phi,n)
   h=list(map(lambda x: phi**(x-1)*(1/x-phi/(x+1)), range(1,n+1)))
   h.insert(0,-1)
   p1=seq_convolute(p,h)
   p2=seq_convolute(p1,h)
   score=[(p1/p)[i] for i in data]
   info=[((p1/p)**2-p2/p)[i] for i in data]
   return [sum(score),sum(info)]

# Jan 28, 2020: Jones median estimator
def jones_median_est (data):
   r=np.median(data)
   est=(r-0.693147)/(np.log(r)+0.3665)
   return est

# Jan 28, 2020: the P0 method
def LD_p0_est (data):
   p0=list(data).count(0)/len(data)
   if (p0>0):
      return -np.log(p0)
   else:
      print("The P0 method is not applicable.")





















