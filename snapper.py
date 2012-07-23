#
#  snapper.py
#  
#  Joint zero-inflated hierarchcial negative binomial model of NZ snapper
#  PyMC implementation of Smith et al. (2012) Ecology: http://www.esajournals.org/doi/abs/10.1890/12-0460.1
# 
#  Created by M. Aaron MacNeil on 20/07/12.
#

from pymc import *
import numpy as np
from scipy.special import gamma as gammaf
from scipy.misc import factorial as ft
from data import *


# = = = = = = = = = = = = Priors
# Overdispersion for NB
delta = Gamma('delta', alpha=0.0001, beta=0.0001, value=6.)
# Zero link intercept to count model
gamma0 = Normal('gamma0', mu=0.0, tau=0.01, value=-0.5)
# Zero link slope to count model
gamma1 = Normal('gamma1', mu=0.0, tau=0.01, value=-0.4)

# Count intercept
beta0 = Normal('beta0', mu=0.0, tau=0.01, value=0.1)
# Reserve effect
beta_reserve = Normal('beta_reserve', mu=0.0, tau=0.01, value=0.)
# Constrain reserve effects
beta_rez = Lambda('beta_rez', lambda b0=beta0,b1=beta_reserve: array([b0-b1,b0+b1]))
# Season effect
beta_season = Normal('beta_season', mu=0.0, tau=0.01, value=0.)
# Area effect
sd_area = Uniform('sd_area', lower=0, upper=5, value=0.4)
tau_area = Lambda('tau_area', lambda sd=sd_area: sd**-2)
beta_area = Normal('beta_area', mu=0.0, tau=tau_area, value=np.zeros(6))
# Reserve given area
beta_rez_area = Lambda('beta_rez_area', lambda br=beta_rez[Ir],ba=beta_area: br+ba)
# Year effect
sd_year = Uniform('sd_year', lower=0, upper=5, value=0.4)
tau_year = Lambda('tau_year', lambda sd=sd_year: sd**-2)
beta_year = Normal('beta_year', mu=0.0, tau=tau_year, value=np.zeros(9))


# = = = = = = = = = = = = Likelihood
# Count model
eta = Lambda('eta', lambda b0=beta_rez_area[Ia],b1=beta_season,b2=beta_year[Iy]: b0+b1*season+b2, trace=False)
# Link function for negative binomial
lambduh = Lambda('lambduh', lambda e=eta: exp(e), trace=False)
# Zero model
pzero = Lambda('pzero', lambda g0=gamma0,g1=gamma1,eta=eta: invlogit(g0+g1*eta), trace=False)

# ZINB likelihood
@observed(dtype=int, plot=False)
def zinb(value=sna, mu=lambduh, alpha=delta, psi=pzero):
    # Initialise likeihood
    like = 0.0
    # Add zero component; zero probability + P(NB==0); value flags for non-zeros to cancel out
    like += sum((np.log(psi + (1.-psi)*(alpha/(mu+alpha))**alpha))*Iz)
    # Add count component; non-zero probability + P(NB>0); value flags for zeros to cancel out
    like += sum((np.log(1.-psi) + np.log(gammaf(alpha+value))-np.log((ft(value)*gammaf(alpha))) + alpha*np.log(alpha/(mu+alpha)) + value*np.log(mu/(mu+alpha)))*Ic)
    return like


# = = = = = = = = = = = = Posteriors
# Zero-inflation probability outside reserve
P_nores = Lambda('P_nores', lambda g0=gamma0,g1=gamma1,bra=beta_rez_area[0]: 1/(1+exp(-(g0+g1*bra))) )
# Zero-inflation probability inside reserve
P_res = Lambda('P_res', lambda g0=gamma0,g1=gamma1,bra=beta_rez_area[1]: 1/(1+exp(-(g0+g1*bra))) )
# Mean count outside reserve
mean_nores = Lambda('mean_nores', lambda pno=P_nores,bra=beta_rez_area[0]: (1-pno)*exp(bra))
# Mean count inside reserve
mean_res = Lambda('mean_res', lambda pno=P_res,bra=beta_rez_area[1]: (1-pno)*exp(bra))
# Reserve effect
res_effect = Lambda('res_effect', lambda no=mean_nores,res=mean_res: res/no)


