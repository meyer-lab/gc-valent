import pymc3 as pm
import theano.tensor as T
import numpy as np
from .model import getTotalActiveSpecies, receptor_expression
from .differencing_op import runCkineDoseOp
from .imports import import_Rexpr, import_pstat


def commonTraf():
    """ Set the common trafficking parameter priors. """
    kfwd = pm.Lognormal('kfwd', mu=np.log(0.001), sd=0.5, shape=1)
    endo = pm.Lognormal('endo', mu=np.log(0.1), sd=0.1, shape=1)
    activeEndo = pm.Lognormal('activeEndo', sd=0.1, shape=1)
    kRec = pm.Lognormal('kRec', mu=np.log(0.1), sd=0.1, shape=1)
    kDeg = pm.Lognormal('kDeg', mu=np.log(0.01), sd=0.2, shape=1)
    sortF = pm.Beta('sortF', alpha=12, beta=80, shape=1)
    return kfwd, endo, activeEndo, kRec, kDeg, sortF

class pSTAT_activity:
    """ Calculating the pSTAT activity residuals for IL2 and IL15 stimulation in distinct cell lines from Visterra. """