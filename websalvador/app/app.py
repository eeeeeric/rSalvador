### Creating a New Post, Dec 21, 2020

# first go, for the classic Luria-Delbruck model

import os

from flask import Flask, render_template, request

import chum as sal

app = Flask(__name__)
app.config['SECRET_KEY'] = os.urandom(24)


@app.route('/')
def index():
    return render_template('index.html', est_m='NA', m_low='NA', m_up='NA',
                           est_mu='NA', mu_low='NA', mu_up='NA', warning=False)


@app.route('/', methods=['POST'])
def my_form_post():
    Nt = request.form['Nt']

    init_m = request.form['init_m']

    Nt = float(Nt)

    init_m = float(init_m)

    mutant0 = request.form['mutantCount']

    mylist = mutant0.split()

    mutants = list(map(lambda y: int(y), mylist))

    est_m, m_low, m_up, est_mu, mu_low, mu_up = [0, 0, 0, 0, 0, 0]

    if init_m < 10:
        warning = False
        est_m = sal.newton_LD(mutants, init_m=init_m)
        m_low, m_up = sal.confint_LD(mutants, init_m=init_m)
        est_mu, mu_low, mu_up = [est_m / Nt, m_low / Nt, m_up / Nt]
        return render_template('indexa.html', est_m=round(est_m, 4),
                               m_low=round(m_low, 4),
                               m_up=round(m_up, 4),
                               est_mu='{:0.2e}'.format(est_mu),
                               mu_low='{:0.2e}'.format(mu_low),
                               mu_up='{:0.2e}'.format(mu_up),
                               rawdata=str(mutant0), given_Nt=str(round(Nt)),
                               given_init_m=str(float(init_m)), warning=warning)

    else:
        warning = True
        return render_template('indexa.html', est_m='NA', m_low='NA',
                               m_up='NA', est_mu='NA', mu_low='NA',
                               mu_up='NA', rawdata=str(mutant0),
                               given_Nt=str(round(Nt)),
                               given_init_m=str(float(init_m)), warning=warning)
