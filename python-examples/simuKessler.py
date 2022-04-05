# Page 792 of Kessler and Levine (2015), eq (48), bw is birth rate of wile-type, etc
# Aug 12, 2021
import numpy as np

def haldane(g=12,N0=30,d=0.1,mu=1e-6):
   Wg=N0
   Mg=0
   for k in range(1,g):
      W_new=np.random.binomial(Wg,1-d)
      M_new=np.random.binomial(W_new, mu)
      Mg=Mg*2+M_new
      Wg=Wg*2-M_new
#      print([Wg, Mg])
   return([Wg, Mg])



# Aug 13, 2021, Markov process version with wild-type cell death
# adapted from simu.Bartlett in rSalvador

def bartlett(b1=1, b2=1, d1=0.2, mu=1e-6, N0=6, T=5.8, maxEvents=1e9, showGrowth=False):
   wild=N0
   mutants=0
   timeNow=0
   events=0
   while (True):
      if showGrowth:
         print("Wild-type cells == "+str(wild)+" and mutant cells == "+str(mutants))
      R=(b1+d1+mu)*wild+b2*mutants
      if R<=0:
         return([-9,-9])
      timeStep=-np.log(np.random.uniform(0,1))/R
      timeNow=timeNow+timeStep
      if (timeNow > T):
         return([wild, mutants])
      else:
         events=events+1
      u_rand=np.random.uniform(0,1)
      if (u_rand <= b1*wild/R):
         wild=wild+1
      elif u_rand<=(b1+d1)*wild/R:
         wild=wild-1
      else:
         mutants=mutants+1 
      if (events > maxEvents):
         print('The number of simulated events has exceeded the limit of '+str(maxEvents))
         return([-9,-9])


# Feb 12, 2022, adding a counter of mutations, extending the bartlett function

def bartlett_mutCounter(b1=1, b2=1, d1=0.2, mu=1e-6, N0=6, T=5.8, maxEvents=1e9, showGrowth=False):
   wild=N0
   mutants=0
   timeNow=0
   events=0
   nMutation=0
   nDivision=0
   while (True):
      if showGrowth:
         print("Wild-type cells == "+str(wild)+" and mutant cells == "+str(mutants))
      R=(b1+d1+mu)*wild+b2*mutants
      if R<=0:
         return([-9,-9, -9, -9])
      timeStep=-np.log(np.random.uniform(0,1))/R
      timeNow=timeNow+timeStep
      if (timeNow > T):
         return([wild, mutants, nDivision, nMutation])
      else:
         events=events+1
      u_rand=np.random.uniform(0,1)
      if (u_rand <= b1*wild/R):
         wild=wild+1
         nDivision+=1
      elif u_rand<=(b1+d1)*wild/R:
         wild=wild-1
      elif u_rand<=(b1+d1+mu)*wild/R:
         mutants=mutants+1
         nDivision+=1
         nMutation=nMutation+1
      else:
         mutants=mutants+1 
      if (events > maxEvents):
         print('The number of simulated events has exceeded the limit of '+str(maxEvents))
         return([-9,-9, -9, -9])


# Aug 16, 2021, change to fixed Nt, modified Aug 18, 2021

def bartlettFixedNt(b1=1, b2=1, d1=0.2, mu=1e-6, N0=6, Nt=1000, maxEvents=1e9, showGrowth=False):
   wild=N0
   mutants=0
   timeNow=0
   events=0
   while (True):
      if showGrowth:
         print("Wild-type cells == "+str(wild)+" and mutant cells == "+str(mutants))
      R=(b1+d1+mu)*wild+b2*mutants
      if R<=0:
         return([-9,-9])
      timeStep=-np.log(np.random.uniform(0,1))/R
      timeNow=timeNow+timeStep
      if (wild >= Nt):
         return([wild, mutants])
      else:
         events=events+1
      u_rand=np.random.uniform(0,1)
      if (u_rand <= b1*wild/R):
         wild=wild+1
      elif u_rand<=(b1+d1)*wild/R:
         wild=wild-1
      else:
         mutants=mutants+1 
      if (events > maxEvents):
         print('The number of simulated events has exceeded the limit of '+str(maxEvents))
         return([-9,-9])

# --------------------------------------

# Aug 19, 2021, adding death of mutants to have a full Bartlett model

def full_bartlettFixedNt(b1=1, b2=1, d1=0.2, d2=0.1, mu=1e-6, N0=6, Nt=1000, maxEvents=1e9, showGrowth=False):
   wild=N0
   mutants=0
   timeNow=0
   events=0
   while (True):
      if showGrowth:
         print("Wild-type cells == "+str(wild)+" and mutant cells == "+str(mutants))
      R=(b1+d1+mu)*wild+(b2+d2)*mutants
      if R<=0:
         print('Population becomes extinct ...')
         return([-9,-9])
      timeStep=-np.log(np.random.uniform(0,1))/R
      timeNow=timeNow+timeStep
      if (wild >= Nt):
         return([wild, mutants])
      else:
         events=events+1
      u_rand=np.random.uniform(0,1)
      if (u_rand <= b1*wild/R):
         wild=wild+1
      elif u_rand<=(b1+d1)*wild/R:
         wild=wild-1
      elif u_rand<=( (b1+d1+mu)*wild + b2*mutants )/R:
         mutants=mutants+1
      else:
         if mutants > 0:
            mutants=mutants-1 
      if (events > maxEvents):
         print('The number of simulated events has exceeded the limit of '+str(maxEvents))
         return([-9,-9])



