"""
This file includes the classes and functions necessary to fit the IL2 and IL15 model to the experimental data.
"""
from os.path import join, dirname, abspath
import pymc3 as pm
import theano.tensor as T
import numpy as np
import pandas as pds
from .model import getTotalActiveSpecies, getSurfaceIL2RbSpecies, getSurfaceGCSpecies
from .differencing_op import runCkineDoseOp


def load_data(filename):
    """ Return path of CSV files. """
    path = dirname(abspath(__file__))
    return pds.read_csv(join(path, filename)).values


class build_model:
    """ Build the overall model handling Ring et al. """

    def __init__(self, traf=True):
        self.traf = traf
        self.M = self.build()

    def build(self):
        """The PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point."""
        # WARNING -- THIS IS JUST THE STUB OF FITTING FOR DATA IMPORT PURPOSES -- WARNING
        M = pm.Model()

        with M:
            endo = pm.Lognormal('endo', mu=np.log(0.1), sd=0.1, shape=1)
            activeEndo = pm.Lognormal('activeEndo', sd=0.1, shape=1)
            kRec = pm.Lognormal('kRec', mu=np.log(0.1), sd=0.1, shape=1)
            kDeg = pm.Lognormal('kDeg', mu=np.log(0.01), sd=0.2, shape=1)
            sortF = pm.Beta('sortF', alpha=12, beta=80, shape=1)
            kfwd = pm.Lognormal('kfwd', mu=np.log(0.001), sd=0.5, shape=1)

            rxnrates = pm.Lognormal('rxn', sd=0.5, shape=6)  # 6 reverse rxn rates for IL2/IL15
            nullRates = T.ones(4, dtype=np.float64)  # k27rev, k31rev, k33rev, k35rev
            Rexpr_2Ra = pm.Lognormal('Rexpr_2Ra', sd=0.5, shape=1)  # Expression: IL2Ra
            Rexpr_2Rb = pm.Lognormal('Rexpr_2Rb', sd=0.5, shape=1)  # Expression: IL2Rb
            Rexpr_15Ra = pm.Lognormal('Rexpr_15Ra', sd=0.5, shape=1)  # Expression: IL15Ra
            Rexpr_gc = pm.Lognormal('Rexpr_gc', sd=0.5, shape=1)  # Expression: gamma chain
            scale = pm.Lognormal('scales', mu=np.log(100.), sd=1, shape=1)  # create scaling constant for activity measurements

            unkVec = T.concatenate((kfwd, rxnrates, nullRates, endo, activeEndo, sortF, kRec, kDeg, Rexpr_2Ra, Rexpr_2Rb, Rexpr_gc, Rexpr_15Ra, nullRates * 0.0))
            unkVec_2Ra_minus = T.concatenate((kfwd, rxnrates, nullRates, endo, activeEndo, sortF, kRec, kDeg, T.zeros(1, dtype=np.float64), Rexpr_2Rb, Rexpr_gc, Rexpr_15Ra, nullRates * 0.0))

            Y_15 = self.dst15.calc(unkVec, scale)  # fitting the data based on dst15.calc for the given parameters
            sd_15 = T.minimum(T.std(Y_15), 0.03)  # Add bounds for the stderr to help force the fitting solution
            pm.Deterministic('Y_15', T.sum(T.square(Y_15)))
            pm.Normal('fitD_15', sd=sd_15, observed=Y_15)  # experimental-derived stderr is used

            Y_int = self.IL2Rb.calc(unkVec)  # fitting the data based on IL2Rb surface data
            sd_int = T.minimum(T.std(Y_int), 0.02)  # Add bounds for the stderr to help force the fitting solution
            pm.Deterministic('Y_int', T.sum(T.square(Y_int)))
            pm.Normal('fitD_int', sd=sd_int, observed=Y_int)

            Y_gc = self.gc.calc(unkVec_2Ra_minus)  # fitting the data using IL2Ra- cells
            sd_gc = T.minimum(T.std(Y_gc), 0.02)  # Add bounds for the stderr to help force the fitting solution
            pm.Deterministic('Y_gc', T.sum(T.square(Y_gc)))
            pm.Normal('fitD_gc', sd=sd_gc, observed=Y_gc)

            # Save likelihood
            pm.Deterministic('logp', M.logpt)

        return M
