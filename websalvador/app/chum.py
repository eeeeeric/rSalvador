### Dec 12, 2020, possibly webSalvador

import rpy2

from rpy2.robjects.packages import importr

base = importr('base')

import rpy2.robjects.packages as rpackages

import rpy2.robjects as robjects

rpackages.importr('rsalvador')

# --------------------------------

demerec_data = list(map(lambda y: int(y), robjects.r['demerec.data']))

foster_data = list(map(lambda y: int(y), robjects.r['foster.data']))

ford_data = list(map(lambda y: int(y), robjects.r['ford.data']))

wh_temp = list(robjects.r['wh.data'])

wh_data = [map(lambda y: int(y), wh_temp[i]) for i in range(0, 13)]


# --------------------------------


def simu_cultures(n=5, mu=1e-8, b1=1, b2=1, N0=5, Nt=1e8):
    ''' simu_cultures(n, mu, b1, b2, N0, Nt) produces n mutant counts. '''

    m = robjects.r['simu.cultures'](n, mu, b1, b2, N0, Nt)

    return (list(map(lambda y: int(y), m)))


def newton_LD(x, phi=1, tol=1e-9, init_m=0, max_iter=30, show_iter=False):
    ''' newton_LD(x, phi=1, tol=1e-9, init_m=0, max_iter=30, show_iter=False) computes
   maximum likelihood estimate of m based on data vector x. '''

    m = robjects.r['newton.LD'](robjects.IntVector(x), phi, tol, init_m,
                                max_iter, show_iter)

    return (m[0])


def newton_LD_plating(x, e, tol=1e-9, init_m=0, max_iter=30, show_iter=False):
    ''' newton_LD(x, e, tol=1e-9, init_m=0, max_iter=30, show_iter=False) computes
   maximum likelihood estimate of m based on data vector x when the plating efficiency
   is e [0 < e < 1.0]. '''

    m = robjects.r['newton.LD.plating'](robjects.IntVector(x), e, tol, init_m,
                                        max_iter, show_iter)

    return (m[0])


# ----------- Mandelbro-Koch model, w is fixed --------------------

def newton_MK(x, w=1, tol=1e-8, init_m=0, max_iter=30, show_iter=False):
    ''' newton_MK(x, w=1, tol=1e-9, init_m=0, max_iter=30, show_iter=False) computes
   maximum likelihood estimate of m based on data vector x when the mutants have a
   relative fitness of w (w>0). It is called a Mandelbrot_Koch model.'''

    m = robjects.r['newton.MK'](robjects.IntVector(x), w, tol, init_m, max_iter,
                                show_iter)

    return (m[0])


def confint_MK(x, w=1, alpha=0.05, tol=1e-8, init_m=0, init_lower=0, init_up=0,
               max_iter=30, show_iter=False):
    '''
   confint_MK(x, w=1, alpha=0.05,tol=1e-9, init_m=0, init_lower=0, init_up=0, max_iter=30, show_iter=False)
   computes a (1-alpha)*100 percet level confidence interval for m, using data vector x, when the mutants
   have a relative fitness of w (w>0).
   '''

    ci = robjects.r['confint.MK'](robjects.IntVector(x), w, alpha, tol, init_m,
                                  init_lower, init_up,
                                  max_iter, show_iter)

    return ([ci[0], ci[1]])


def confint_LD(x, alpha=0.05, phi=1, tol=1e-9, init_m=0, init_lower=0,
               init_up=0, max_iter=30, show_iter=False):
    '''
   confint_LD(x, alpha=0.05, phi=1, tol=1e-9, init_m=0, init_lower=0, init_up=0, max_iter=30, show_iter=False)
   computes a (1-alpha)*100 percet level confidence interval for m, using data vector x.
   '''
    ci = robjects.r['confint.LD'](robjects.IntVector(x), alpha, phi, tol,
                                  init_m, init_lower, init_up,
                                  max_iter, show_iter)

    return ([ci[0], ci[1]])


def confint_LD_plating(x, e, alpha=0.05, tol=1e-9, init_m=0, init_lower=0,
                       init_up=0, max_iter=30, show_iter=False):
    '''
   confint_LD_plating(x, e, alpha=0.05,tol=1e-9, init_m=0, init_lower=0, init_up=0, max_iter=30, show_iter=False)
   computes a (1-alpha)*100 percet level confidence interval for m, using data vector x, when the plating
   efficiency is e [0 < e < 1.0]
   '''

    ci = robjects.r['confint.LD.plating'](robjects.IntVector(x), e, alpha, tol,
                                          init_m, init_lower, init_up,
                                          max_iter, show_iter)

    return ([ci[0], ci[1]])


# ------- likelihood ratio test --------
def LRT_LD(x1, x2, R=1, phi1=1, phi2=1, init_mc=0, init_m1=0, init_m2=0,
           tol=1e-9, max_iter=30, show_iter=False):
    '''
   LRT_LD(x1, x2, R=1, phi1=1, phi2=1, init_mc=0, init_m1=0, init_m2=0, tol=1e-9, max_iter=30, show_iter=False) computes
   a likelihood ratio statistic and the corresponding p-value, for comparing the mutation rates in two experiment. Note
   that R=(Nt in the second experiment) / (Nt in the first experimetn).
   '''

    out = robjects.r['LRT.LD'](robjects.IntVector(x1), robjects.IntVector(x2),
                               R, phi1, phi2, init_mc, init_m1,
                               init_m2, tol, max_iter, show_iter)

    return ([out[0], out[1]])


def LRT_LD_plating(x1, x2, R=1, e1=0.1, e2=0.1, init_mc=0, init_m1=0, init_m2=0,
                   tol=1e-9, max_iter=30, show_iter=False):
    '''
   LRT_LD_plating(x1, x2, R=1, e1=1, e2=1, init_mc=0, init_m1=0, init_m2=0, tol=1e-9, max_iter=30, show_iter=False) computes
   a likelihood ratio statistic and the corresponding p-value, for comparing the mutation rates in two experiment. Note
   that R=(Nt in the second experiment) / (Nt in the first experimetn). The plating efficiencies in the two exepriments
   are e1 and e2, respetively.
   '''

    out = robjects.r['LRT.LD.plating'](robjects.IntVector(x1),
                                       robjects.IntVector(x2), R, e1, e2,
                                       init_mc, init_m1,
                                       init_m2, tol, max_iter, show_iter)

    return ([out[0], out[1]])


def confint_foldchange_MK(x, y, Nx, Ny, w1=1, w2=1, alpha=0.05, init_base=1e-8,
                          init_fold=1.6,
                          init_low_base=-9, init_up_base=-9, init_low_fold=-9,
                          init_up_fold=-9,
                          max_iter=30, tol=1e-9, show_iter=False):
    ci = robjects.r['confint.foldchange.MK'](robjects.IntVector(x),
                                             robjects.IntVector(y), Nx, Ny, w1,
                                             w2,
                                             alpha, init_base, init_fold,
                                             init_low_base, init_up_base,
                                             init_low_fold,
                                             init_up_fold, max_iter, tol,
                                             show_iter)
    return ([ci[0], ci[1], ci[2]])


def boot_foldchange_MK(x, y, N1, N2, w1=1, w2=1, alpha=0.05, init_m1=0,
                       init_m2=0, na_ok=0, nboot=10):
    ci = robjects.r['boot.foldchange.MK'](robjects.IntVector(x),
                                          robjects.IntVector(y), N1, N2, w1, w2,
                                          alpha, init_m1, init_m2, na_ok, nboot)
    #  return([ci[0], ci[1],ci[2]])
    return (ci)


def bayes_foldchange_MK(x, y, N1, N2, w1=1, w2=1, s0=0.4, s1=0.5,
                        init_base=1e-8, init_fold=1.6,
                        v0=100, v1=100, iter=1100, burn=100, thin=1, alpha=0.05,
                        short_out=True, show_simu=False):
    ci = robjects.r['bayes.foldchange.MK'](robjects.IntVector(x),
                                           robjects.IntVector(y), N1, N2, w1,
                                           w2,
                                           s0, s1, init_base, init_fold, v0, v1,
                                           iter, burn, thin, alpha, short_out,
                                           show_simu)
    return ([ci[0], ci[1], ci[2]])
