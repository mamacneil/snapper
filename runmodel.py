# Run model
#
# PyMC implementation of Smith et al. (2012) Ecology: http://www.esajournals.org/doi/abs/10.1890/12-0460.1
# 
#  Created by M. Aaron MacNeil on 20/07/12.
#

import snapper
from pymc import MCMC, BinaryMetropolis, Metropolis, AdaptiveMetropolis
from pymc import Matplot as mp


M = MCMC(models)
xex = 6
M.isample(10**xex, 10**xex-10**(xex-1), thin=100, verbose=2)

try:
    os.mkdir('Outputs')
except OSError:
    pass
os.chdir('Outputs')
M.write_csv("results.csv")
mp.plot(M)
os.chdir('..')