import mpmath as mp
mp.mp.dps = 40

# Apr 7, 2022, plating+fitness, Kessler-Levine approach

def kessler_integrand(t,k,m,w,e):
   theta=(1-e)/e
   a1=m*mp.gamma(1/w -1)/w/mp.gamma(1/w) *t/(t+theta+1)\
   *mp.hyp2f1(1,1-1/w,2-1/w, t/(t+theta+1))
   a2=-m*mp.pi/w*(t/(t+theta+1))**(1/w)
   a=mp.exp(a1+a2/mp.tan(mp.pi/w))*mp.sin(-a2)/(1+t)**(k+1)
   return a


def prob_kessler(k=6,m=1.3,w=1.2, e=0.1):
   prob= mp.quad(lambda t:kessler_integrand(t, k=k, m=m, w=w, e=e), [0, mp.inf])
   return(float(prob/mp.pi))

# Apr 8, 2022

def sharp_p0(m, w, e):
   theta=(1-e)/e
   p=-m*mp.hyp2f1(1,1/w,1+1/w, -theta)
   return float(mp.exp(p))

def sharp_p1(m, w, e):
   theta=(1-e)/e
   p0=sharp_p0(m,w,e)
   p=m*p0*mp.hyp2f1(2,1+1/w,2+1/w, -theta)/e/(1+w)
   return float(p)

# Apr 9, 2022, exact p2

def sharp_p2(m, w, e):
   theta=(1-e)/e
   p0=sharp_p0(m,w,e)
   h1=mp.hyp2f1(2,1+1/w,2+1/w, -theta)
   h2=mp.hyp2f1(3,2+1/w,3+1/w, -theta)
   f1=(1+theta)**2 *m**2 *h1**2 /(1+w)**2
   f2=2*m*(1+theta)*h1/(1+w)
   f3=2*m*(1+theta)**2 * h2/(1+2*w)
   p=0.5*p0*(f1+f2-f3)
   return float(p)
