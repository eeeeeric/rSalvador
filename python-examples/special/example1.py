#Sep 1, 2022, demo for the fresh approach paper

import numpy as np

import gasch_likely as zheng

import gasch_confint as daisy

from datetime import datetime
now0=datetime.now()


Nt=[881200,1147200, 529800,1215300, 230000, 748400, 296500, 378800,1318500, 1328000, 999400,1567500]


plat=[0.12, 0.11, 0.22, 0.14, 0.20, 0.04, 0.40, 0.87, 0.63, 0.27, 0.28, 0.50]


ee=np.array(plat)/100

ww=[0.8]*12

mutant=[2, 1, 19, 42, 10,  0,  6,  8, 32, 10,  3, 11]


point=zheng.golden_plafit(data=mutant, Nt=Nt, w=ww, e=ee, mu_low=1e-6, mu_up=0.01)

ci=daisy.confint_bisect(data=mutant, Nt=Nt, w=ww, e=ee, mu_low=1e-6, mu_up=0.01)

now=datetime.now()

print('Work ending at ...',now)
print('The total amount of time consumed is ...',now-now0, ' goodbye')



