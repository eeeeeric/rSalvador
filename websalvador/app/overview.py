### Creating a New Post, Dec 21, 2020

# first go, for the classic Luria-Delbruck model

import sys, os

from flask import Flask, render_template, request

import numpy as np

import chum as sal

from useful import prune, is_number

LARGE_COUNT=8000

app=Flask(__name__)
app.config['SECRET_KEY'] = os.urandom(24)


@app.route('/')

def index():

   return render_template('dome.html')


## ----------- Basic Lee-Coulson -----------------


@app.route('/lc_basic')

def index_lc_basic():

   return render_template('loquat6.html', est_m='NA', m_low='NA', m_up='NA',\
         est_mu='NA', mu_low='NA', mu_up='NA', warning=False)

@app.route('/lc_basic', methods=['POST'])
def my_form_post():

   Nt=request.form['Nt']
   if not is_number(Nt):
      Nt=1.9e8
   init_m=request.form['init_m']
   if not is_number(init_m):
      init_m=1.58
   init_m_low=request.form['init_m_low']
   if not is_number(init_m_low):
      init_m_low=0.58
   init_m_up=request.form['init_m_up']
   if not is_number(init_m_up):
      init_m_up=2.58

   Nt=float(Nt)
   init_m=float(init_m)
   init_m_low=float(init_m_low)
   init_m_up=float(init_m_up)

   mutant0 = request.form['mutantCount']

   if not all(list(map(lambda y:y.isdigit(), mutant0.split()))):

       warning=True
       return render_template('loquat6a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
      
   mylist=mutant0.split()

   mutants=list(map(lambda y:int(y), mylist))

   mutants=list(map(lambda z:prune(z,LARGE_COUNT),mutants))

   try:
       est_m=sal.newton_LD(mutants, init_m=init_m)
       m_low, m_up =sal.confint_LD(mutants, init_m=init_m, init_lower=init_m_low, init_up=init_m_up)
       est_mu, mu_low, mu_up =[est_m/Nt, m_low/Nt, m_up/Nt]

   except:
       warning=True
       return render_template('loquat6a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
   else:
      warning=False   
      return render_template('loquat6a.html', est_m=round(est_m,4), m_low=round(m_low,4),\
          m_up=round(m_up,4), est_mu='{:0.2e}'.format(est_mu), mu_low='{:0.2e}'.format(mu_low),\
          mu_up='{:0.2e}'.format(mu_up), rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )




## ----------- Lee-Coulson with an incomplating plating efficiency -----------------


@app.route('/lc_plating')

def index_lc_plating():

   return render_template('mulbery1.html', est_m='NA', m_low='NA', m_up='NA',\
         est_mu='NA', mu_low='NA', mu_up='NA', warning=False)

@app.route('/lc_plating', methods=['POST'])
def lc_platin_post():

   Nt=request.form['Nt']
   if not is_number(Nt):
      Nt=1.9e8
   e=request.form['e']
   if not is_number(e):
      e=0.1
   if (float(e)<=0) or (float(e)>=1):
      e=0.1
   init_m=request.form['init_m']
   if not is_number(init_m):
      init_m=1.58
   init_m_low=request.form['init_m_low']
   if not is_number(init_m_low):
      init_m_low=0.58
   init_m_up=request.form['init_m_up']
   if not is_number(init_m_up):
      init_m_up=2.58



   Nt=float(Nt)
   e=float(e)
   init_m=float(init_m)
   init_m_low=float(init_m_low)
   init_m_up=float(init_m_up)

   mutant0 = request.form['mutantCount']

   if not all(list(map(lambda y:y.isdigit(), mutant0.split()))):

       warning=True
       return render_template('mulbery1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
      
   mylist=mutant0.split()

   mutants=list(map(lambda y:int(y), mylist))
   mutants=list(map(lambda z:prune(z,LARGE_COUNT),mutants))

   try:
       est_m=sal.newton_LD_plating(mutants, e=e,init_m=init_m)
       m_low, m_up =sal.confint_LD_plating(mutants, e=e, init_m=init_m, init_lower=init_m_low, init_up=init_m_up)
       est_mu, mu_low, mu_up =[est_m/Nt, m_low/Nt, m_up/Nt]

   except:
       warning=True
       return render_template('mulbery1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_e=str(float(e)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
   else:
      warning=False   
      return render_template('mulbery1a.html', est_m=round(est_m,4), m_low=round(m_low,4),\
          m_up=round(m_up,4), est_mu='{:0.2e}'.format(est_mu), mu_low='{:0.2e}'.format(mu_low),\
          mu_up='{:0.2e}'.format(mu_up), rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_e=str(float(e)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )




## ----------- Lee-Coulson with a varying Nt -----------------


@app.route('/lc_varyingNt')

def index_lc_Nt():

   return render_template('egret1.html', est_m='NA', m_low='NA', m_up='NA',\
         est_mu='NA', mu_low='NA', mu_up='NA', warning=False)

@app.route('/lc_varyingNt', methods=['POST'])
def lc_varyingNt_post():
   
   Nt=request.form['Nt']
   if not is_number(Nt):
      Nt=1.9e8
   cv=request.form['cv']
   if not is_number(cv):
      cv=0.2
   if (float(cv)<=0) or (float(cv)>=5):
      cv=0.2
   init_m=request.form['init_m']
   if not is_number(init_m):
      init_m=1.58
   init_m_low=request.form['init_m_low']
   if not is_number(init_m_low):
      init_m_low=0.58
   init_m_up=request.form['init_m_up']
   if not is_number(init_m_up):
      init_m_up=2.58


   Nt=float(Nt)
   cv=float(cv)
   init_m=float(init_m)
   init_m_low=float(init_m_low)
   init_m_up=float(init_m_up)

   mutant0 = request.form['mutantCount']

   if not all(list(map(lambda y:y.isdigit(), mutant0.split()))):

       warning=True
       return render_template('egret1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
      
   mylist=mutant0.split()

   mutants=list(map(lambda y:int(y), mylist))
   mutants=list(map(lambda z:prune(z,LARGE_COUNT),mutants))

   try:
       est_m=sal.newton_B0(mutants, cv=cv,init_m=init_m)
       m_low, m_up =sal.confint_B0(mutants, cv=cv, init_m=init_m, init_lower=init_m_low, init_up=init_m_up)
       est_mu, mu_low, mu_up =[est_m/Nt, m_low/Nt, m_up/Nt]

   except:
       warning=True
       return render_template('egret1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_cv=str(float(cv)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
   else:
      warning=False   
      return render_template('egret1a.html', est_m=round(est_m,4), m_low=round(m_low,4),\
          m_up=round(m_up,4), est_mu='{:0.2e}'.format(est_mu), mu_low='{:0.2e}'.format(mu_low),\
          mu_up='{:0.2e}'.format(mu_up), rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_cv=str(float(cv)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )



## ----------- Mandelbrot-Koch model, Jan 3, 2020 -------------


@app.route('/mk')

def index_mk():

   return render_template('cherry1.html', est_m='NA', m_low='NA', m_up='NA',\
         est_mu='NA', mu_low='NA', mu_up='NA', warning=False)

@app.route('/mk', methods=['POST'])
def mk_form_post():
   
   Nt=request.form['Nt']
   if not is_number(Nt):
      Nt=1.9e8
   w=request.form['w']
   if not is_number(w):
      w=1.0
   if (float(w)<=0) or (float(w)>=3):
      w=1.0
   init_m=request.form['init_m']
   if not is_number(init_m):
      init_m=1.58
   init_m_low=request.form['init_m_low']
   if not is_number(init_m_low):
      init_m_low=0.58
   init_m_up=request.form['init_m_up']
   if not is_number(init_m_up):
      init_m_up=2.58


   Nt=float(Nt)
   w=float(w)
   init_m=float(init_m)
   init_m_low=float(init_m_low)
   init_m_up=float(init_m_up)

   mutant0 = request.form['mutantCount']

   if not all(list(map(lambda y:y.isdigit(), mutant0.split()))):

       warning=True
       return render_template('cherry1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
      
   mylist=mutant0.split()

   mutants=list(map(lambda y:int(y), mylist))
   mutants=list(map(lambda z:prune(z,LARGE_COUNT),mutants))

   try:
       est_m=sal.newton_MK(mutants, w=w,init_m=init_m)
       m_low, m_up =sal.confint_MK(mutants, w=w, init_m=init_m, init_lower=init_m_low, init_up=init_m_up)
       est_mu, mu_low, mu_up =[est_m/Nt, m_low/Nt, m_up/Nt]

   except:
       warning=True
       return render_template('cherry1a.html', est_m='NA', m_low='NA',\
          m_up='NA', est_mu='NA', mu_low='NA',\
          mu_up='NA', rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_w=str(float(w)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )
   else:
      warning=False   
      return render_template('cherry1a.html', est_m=round(est_m,4), m_low=round(m_low,4),\
          m_up=round(m_up,4), est_mu='{:0.2e}'.format(est_mu), mu_low='{:0.2e}'.format(mu_low),\
          mu_up='{:0.2e}'.format(mu_up), rawdata=str(mutant0), given_Nt=str(round(Nt)),\
          given_w=str(float(w)),\
          given_init_m=str(float(init_m)), given_init_m_low=str(float(init_m_low)),\
          given_init_m_up=str(float(init_m_up)), warning=warning  )








## ----------- comparison under Mandelbrot-Koch, Jan 2, 2020 -----------------


@app.route('/compare_mk')

def index_compare_mk():

   return render_template('kiwi1.html', lrt_stat='NA', p_value='NA', init_m1='NA',\
         init_m2='NA', init_mc='NA',  warning=False)

@app.route('/compare_mk', methods=['POST'])
def compare_mk_form_post():

   Nt1=request.form['Nt_1']
   if not is_number(Nt1):
      Nt1=1.9e8

   Nt2=request.form['Nt_2']
   if not is_number(Nt2):
      Nt2=1.9e8

   w1=request.form['w_1']
   if not is_number(w1):
      w1=1.0

   w2=request.form['w_2']
   if not is_number(w2):
      w2=1.0


   init_m1=request.form['init_m1']
   if not is_number(init_m1):
      init_m1=1.58

   init_m2=request.form['init_m2']
   if not is_number(init_m2):
      init_m2=1.58

   init_mc=request.form['init_mc']
   if not is_number(init_mc):
      init_mc=1.58

   Nt1=float(Nt1)
   Nt2=float(Nt2)
   R=Nt2/Nt1
   w1=float(w1)
   w2=float(w2)

   init_m1=float(init_m1)
   init_m2=float(init_m2)
   init_mc=float(init_mc)

   mutant1 = request.form['mutantCount1']
   mutant2 = request.form['mutantCount2']

   if ( not all(list(map(lambda y:y.isdigit(), mutant1.split())))  ) or\
    (  not  all(list(map(lambda y:y.isdigit(), mutant2.split())))  ):
    
      warning=True
      return render_template('kiwi1a.html', lrt_stat='NA', p_value='NA',\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )

   mylist0=mutant1.split()
   x1=list(map(lambda y:int(y), mylist0))

   mylist0=mutant2.split()
   x2=list(map(lambda y:int(y), mylist0))

   try:

      lrt_stat, p_value =sal.LRT_MK(x1=x1, x2=x2, R=R, w1=w1, w2=w2, init_mc=init_mc,init_m1=init_m1,init_m2=init_m2)

   except:
      warning=True
      return render_template('kiwi1a.html', lrt_stat='NA', p_value='NA',\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )
   else:
      warning=False
      return render_template('kiwi1a.html', lrt_stat=round(lrt_stat,4), p_value='{:0.4e}'.format(p_value),\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )



## ----------- comparison under Lea-Coulson with partial plating, Jan 2, 2020 -----------------


@app.route('/compare_lc_plating')

def index_compare_lc_plating():

   return render_template('mango1.html', lrt_stat='NA', p_value='NA', init_m1='NA',\
         init_m2='NA', init_mc='NA',  warning=False)

@app.route('/compare_lc_plating', methods=['POST'])
def compare_lc_plating_form_post():

   Nt1=request.form['Nt_1']
   if not is_number(Nt1):
      Nt1=1.9e8

   Nt2=request.form['Nt_2']
   if not is_number(Nt2):
      Nt2=1.9e8

   e1=request.form['e_1']
   if not is_number(e1):
      e1=0.1

   e2=request.form['e_2']
   if not is_number(e2):
      e2=0.1


   init_m1=request.form['init_m1']
   if not is_number(init_m1):
      init_m1=1.58

   init_m2=request.form['init_m2']
   if not is_number(init_m2):
      init_m2=1.58

   init_mc=request.form['init_mc']
   if not is_number(init_mc):
      init_mc=1.58


   Nt1=float(Nt1)
   Nt2=float(Nt2)
   R=Nt2/Nt1
   e1=float(e1)
   e2=float(e2)

   init_m1=float(init_m1)
   init_m2=float(init_m2)
   init_mc=float(init_mc)

   mutant1 = request.form['mutantCount1']
   mutant2 = request.form['mutantCount2']

   if ( not all(list(map(lambda y:y.isdigit(), mutant1.split())))  ) or\
    (  not  all(list(map(lambda y:y.isdigit(), mutant2.split())))  ):
    
      warning=True
      return render_template('mango1a.html', lrt_stat='NA', p_value='NA',\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )

   mylist0=mutant1.split()
   x1=list(map(lambda y:int(y), mylist0))

   mylist0=mutant2.split()
   x2=list(map(lambda y:int(y), mylist0))

   try:

      lrt_stat, p_value =sal.LRT_LD_plating(x1=x1, x2=x2, R=R, e1=e1, e2=e2,
               init_mc=init_mc,init_m1=init_m1,init_m2=init_m2)

   except:
      warning=True
      return render_template('mango1a.html', lrt_stat='NA', p_value='NA',\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )
   else:
      warning=False
      return render_template('mango1a.html', lrt_stat=round(lrt_stat,4), p_value='{:0.4e}'.format(p_value),\
          given_init_m1=str(round(init_m1,2)), given_init_m2=str(round(init_m2,2)),\
          given_init_mc=str(round(init_mc,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )




## ----------- fold change under Mandelbrot-Koch, Jan 3, 2020 -----------------

@app.route('/foldchange_mk')

def index_fold_mk():

   return render_template('apricot1.html', fold_est='NA', fold_ci_low='NA', fold_ci_up='NA',\
         init_m2='NA', init_mc='NA',  warning=False)

@app.route('/foldchange_mk', methods=['POST'])
def fold_mk_form_post():

   Nt1=request.form['Nt_1']
   if not is_number(Nt1):
      Nt1=1.9e8

   Nt2=request.form['Nt_2']
   if not is_number(Nt2):
      Nt2=1.9e8

   w1=request.form['w_1']
   if not is_number(w1):
      w1=1.0

   w2=request.form['w_2']
   if not is_number(w2):
      w2=1.0


   init_base=request.form['init_base']
   if not is_number(init_base):
      init_base=1e-7

   init_fold=request.form['init_fold']
   if not is_number(init_fold):
      init_fold=1.0

   init_low_base=request.form['init_low_base']
   if not is_number(init_low_base):
      init_low_base=-9

   init_low_fold=request.form['init_low_fold']
   if not is_number(init_low_fold):
      init_low_fold=-9


   init_up_base=request.form['init_up_base']
   if not is_number(init_up_base):
      init_up_base=-9

   init_up_fold=request.form['init_up_fold']
   if not is_number(init_up_fold):
      init_up_fold=-9



   Nt1=float(Nt1)
   Nt2=float(Nt2)
   w1=float(w1)
   w2=float(w2)

   init_base=float(init_base)
   init_fold=float(init_fold)
   init_low_base=float(init_low_base)
   init_low_fold=float(init_low_fold)
   init_up_base=float(init_up_base)
   init_up_fold=float(init_up_fold)


   mutant1 = request.form['mutantCount1']
   mutant2 = request.form['mutantCount2']

   if ( not all(list(map(lambda y:y.isdigit(), mutant1.split())))  ) or\
    (  not  all(list(map(lambda y:y.isdigit(), mutant2.split())))  ):
    
      warning=True
      return render_template('apricot1a.html', fold_est='NA', fold_ci_low='NA',\
          fold_ci_up='NA',\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )

   mylist0=mutant1.split()
   x1=list(map(lambda y:int(y), mylist0))

   mylist0=mutant2.split()
   x2=list(map(lambda y:int(y), mylist0))

   try:

       fold_ci_low, fold_est, fold_ci_up =sal.confint_foldchange_MK(x=x1, y=x2, Nx=Nt1, Ny=Nt2,\
          w1=w1, w2=w2, init_base=init_base, init_fold=init_fold, init_low_base=init_low_base, \
          init_up_base=init_up_base,\
          init_low_fold=init_low_fold, init_up_fold=init_up_fold)

   except:
      warning=True
      return render_template('apricot1a.html', fold_est='NA', fold_ci_low='NA',\
          fold_ci_up='NA',\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )

   else:
      warning=False
      return render_template('apricot1a.html', fold_est=round(fold_est,4), fold_ci_low=round(fold_ci_low,4),\
          fold_ci_up=round(fold_ci_up,4),\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_w_1=str(round(w1,2)),\
          given_w_2=str(round(w2,2)),warning=warning  )



## ----------- fold change under Mandelbrot-Koch, Jan 3, 2020 -----------------

@app.route('/fold_ld_plating')

def index_fold_ld_plating():

   return render_template('plum1.html', fold_est='NA', fold_ci_low='NA', fold_ci_up='NA',\
         init_m2='NA', init_mc='NA',  warning=False)

@app.route('/fold_ld_plating', methods=['POST'])
def fold_ld_plating_form_post():

   Nt1=request.form['Nt_1']
   if not is_number(Nt1):
      Nt1=1.9e8

   Nt2=request.form['Nt_2']
   if not is_number(Nt2):
      Nt2=1.9e8

   e1=request.form['e_1']
   if not is_number(e1):
      e1=0.1

   e2=request.form['e_2']
   if not is_number(e2):
      e2=0.1


   init_base=request.form['init_base']
   if not is_number(init_base):
      init_base=1e-7

   init_fold=request.form['init_fold']
   if not is_number(init_fold):
      init_fold=1.0

   init_low_base=request.form['init_low_base']
   if not is_number(init_low_base):
      init_low_base=-9

   init_low_fold=request.form['init_low_fold']
   if not is_number(init_low_fold):
      init_low_fold=-9


   init_up_base=request.form['init_up_base']
   if not is_number(init_up_base):
      init_up_base=-9

   init_up_fold=request.form['init_up_fold']
   if not is_number(init_up_fold):
      init_up_fold=-9


   Nt1=float(Nt1)
   Nt2=float(Nt2)
   e1=float(e1)
   if e1<0 or e1>=1:
      e1=0.1
   e2=float(e2)
   if e2<0 or e2>=1:
      e2=0.1

   init_base=float(init_base)
   init_fold=float(init_fold)
   init_low_base=float(init_low_base)
   init_low_fold=float(init_low_fold)
   init_up_base=float(init_up_base)
   init_up_fold=float(init_up_fold)


   mutant1 = request.form['mutantCount1']
   mutant2 = request.form['mutantCount2']

   if ( not all(list(map(lambda y:y.isdigit(), mutant1.split())))  ) or\
    (  not  all(list(map(lambda y:y.isdigit(), mutant2.split())))  ):
    
      warning=True
      return render_template('plum1a.html', fold_est='NA', fold_ci_low='NA',\
          fold_ci_up='NA',\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )

   mylist0=mutant1.split()
   x1=list(map(lambda y:int(y), mylist0))

   mylist0=mutant2.split()
   x2=list(map(lambda y:int(y), mylist0))

   try:

       fold_ci_low, fold_est, fold_ci_up =sal.confint_foldchange_LD_plating(x=x1, y=x2, Nx=Nt1, Ny=Nt2,\
          e1=e1, e2=e2, init_base=init_base, init_fold=init_fold, init_low_base=init_low_base, \
          init_up_base=init_up_base,\
          init_low_fold=init_low_fold, init_up_fold=init_up_fold)

   except:
      warning=True
      return render_template('plum1a.html', fold_est='NA', fold_ci_low='NA',\
          fold_ci_up='NA',\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )

   else:
      warning=False
      return render_template('plum1a.html', fold_est=round(fold_est,4), fold_ci_low=round(fold_ci_low,4),\
          fold_ci_up=round(fold_ci_up,4),\
          given_init_base=str(round(init_base,9)), given_init_fold=str(round(init_fold,2)),\
          given_init_low_base=str(round(init_low_base,2)),\
          given_init_low_fold=str(round(init_low_fold,2)), given_init_up_base=str(round(init_up_base,2)),\
          given_init_up_fold=str(round(init_up_fold,2)),\
          given_x1=str(mutant1), given_x2=str(mutant2),given_Nt_1=str(round(Nt1)),\
          given_Nt_2=str(round(Nt2)), given_e_1=str(round(e1,2)),\
          given_e_2=str(round(e2,2)),warning=warning  )









