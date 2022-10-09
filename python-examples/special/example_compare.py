### Sep 4, 2022, LR test for fresh approach paper

import numpy as np
import sys
import os

import gasch_testing as zheng

from datetime import datetime
now0=datetime.now()

print('Work starting at ...',now0)

# --------------- expt A ------------------
Nt1=[432900, 54300,145600,103700,138600,115000,100100, 51400,364100,11880]

plat=[0.86, 5.61, 2.40, 4.70, 3.69, 5.25, 3.57, 8.14, 1.46, 3.93]

e1=np.array(plat)/100

w1=[1.5]*10

mutant1=[213,  31,  481,   79,  151,  161,  833,  895, 1262,  899]

# --------------- expt B ------------------
Nt2=[881200,1147200, 529800,1215300, 230000, 748400, 296500, 378800,1318500, 1328000, 999400,1567500]

plat=[0.12, 0.11, 0.22, 0.14, 0.20, 0.04, 0.40, 0.87, 0.63, 0.27, 0.28, 0.50]

e2=np.array(plat)/100

w2=[0.8]*12

mutant2=[2, 1, 19, 42, 10,  0,  6,  8, 32, 10,  3, 11]

out=zheng.plafit_LR(mutant1, mutant2, Nt1, Nt2, w1, w2, e1, e2, mu_up1=5e-3, mu_up2=5e-3, mu_upc=5e-3)

now=datetime.now()

print('Work ending at ...',now)
print('The total amount of time consumed is ...',now-now0, ' goodbye')







