#Sep 1, 2022, demo for the fresh paper paper

import numpy as np

import gasch_likely as zheng

import gasch_confint as daisy

from datetime import datetime
now0=datetime.now()


Nt=[432900, 54300,145600,103700,138600,115000,100100, 51400,364100,11880]

plat=[0.86, 5.61, 2.40, 4.70, 3.69, 5.25, 3.57, 8.14, 1.46, 3.93]

ee=np.array(plat)/100

ww=[1.5]*10

mutant=[213,  31,  481,   79,  151,  161,  833,  895, 1262,  899]


#point=zheng.golden_plafit(data=mutant, Nt=Nt, w=ww, e=ee, mu_low=1e-6, mu_up=0.01)

ci=daisy.confint_bisect(data=mutant, Nt=Nt, w=ww, e=ee, mu_low=1e-6, mu_up=0.01)

now=datetime.now()

print('Work ending at ...',now)
print('The total amount of time consumed is ...',now-now0, ' goodbye')



