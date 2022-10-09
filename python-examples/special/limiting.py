## Apr 18, 2022, asymptotic behavior of mutant probability

import math

def limit(k=200, m=2.2, e=0.01, w=0.9):
   asymp=m*math.gamma(1+1/w) * e**(1/w) /k**(1+1/w)  # modified Apr 26, 2022, if w>1.0
   return asymp/w
