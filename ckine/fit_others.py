"""
This file includes the classes and functions necessary to fit the IL4, IL7 and IL21 model to the experimental data.
"""
import pymc3 as pm, theano.tensor as T, os
from os.path import join
from theano import shared
import numpy as np, pandas as pds
from .model import getActiveSpecies, getSurfaceIL2RbSpecies
from .differencing_op import runCkineOp, runCkineKineticOp

class IL4_7_21_activity:
    def __init__(self):
        """This loads the experiment data and saves it as a member matrix and it also makes a vector of the ligand concentrations that we are going to take care of."""
        path = os.path.dirname(os.path.abspath(__file__))
        dataIL4 = pds.read_csv(join(path, "./data/Gonnard-Fig1B.csv")).as_matrix() # imports IL4 file into pandas array
        dataIL21 = pds.read_csv(join(path, "./data/Gonnard-Fig1C.csv")).as_matrix() # imports IL21 file into pandas array
        dataIL7 = pds.read_csv(join(path, "./data/Gonnard-FigS4B.csv")).as_matrix() # imports IL7 file into pandas array
        self.cytokC = np.array([5., 50., 500., 5000., 50000.])

        self.fit_data = np.concatenate((dataIL4[:, 1], dataIL7[:, 1], dataIL21[:, 1])) # the measurements are normalized to 1

        npactivity = getActiveSpecies().astype(np.float64)
        self.activity = shared(np.concatenate((npactivity, 0.5*npactivity, np.zeros(6)))) # 0.5 is because its the endosome

    def calc(self, unkVec):
        """Simulate the experiment with different ligand stimulations. It is making a list of promises which will be calculated and returned as output."""

        Op = runCkineOp(ts=np.array(500.))

        # Loop over concentrations of IL4
        actVecIL4 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[4], x))), self.cytokC)))

        # Loop over concentrations of IL7
        actVecIL7 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[2], x))), self.cytokC)))

        # Loop over concentrations of IL21
        actVecIL21 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[5], x))), self.cytokC)))

        # Normalize to the maximal activity, put together into one vector
        actVec = T.concatenate((actVecIL4 / T.max(actVecIL4), actVecIL7 / T.max(actVecIL7), actVecIL21 / T.max(actVecIL21)))

        # value we're trying to minimize is the distance between the y-values on points of the graph that correspond to the same lignad values and species
        return self.fit_data - actVec


class build_model:
    """Going to load the data from the CSV file at the very beginning of when build_model is called... needs to be separate member function to avoid uploading file thousands of times."""
    def __init__(self):
        self.act = IL4_7_21_activity()
        self.M = self.build()

    def build(self):
        """The PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point."""
        M = pm.Model()

        with M:
            kfwd = pm.Lognormal('kfwd', mu=np.log(0.00001), sd=0.1)
            rxnrates = pm.Lognormal('rxn', mu=np.log(0.1), sd=0.1, shape=10) # size of all unknown off rates in model
            endo_activeEndo = pm.Lognormal('endo', mu=np.log(0.1), sd=0.1, shape=2)
            kRec_kDeg = pm.Lognormal('kRec_kDeg', mu=np.log(0.1), sd=0.1, shape=2)
            GCexpr = pm.Lognormal('GCexpr', sd=0.1, shape=1) # expression of gc
            IL7Raexpr = pm.Lognormal('IL7Raexpr', sd=0.1, shape=1) # expression of IL7Ra
            Rexpr = pm.Lognormal('Rexpr', sd=0.1, shape=2) # Expression: IL4Ra, IL21Ra
            sortF = pm.Beta('sortF', alpha=20, beta=40, testval=0.333)*0.95

            ligands = T.zeros(6, dtype=np.float64)

            unkVec = T.concatenate((ligands, T.stack(kfwd), rxnrates, endo_activeEndo, T.stack(sortF), kRec_kDeg, T.zeros(2, dtype=np.float64), T.stack(GCexpr), T.zeros(1, dtype=np.float64), T.stack(IL7Raexpr), T.zeros(1, dtype=np.float64), T.stack(Rexpr))) # receptor expression indexing same as in model.cpp

            Y_int = self.act.calc(unkVec) # fitting the data based on act.calc for the given parameters

            pm.Deterministic('Y_int', T.sum(T.square(Y_int)))

            pm.Normal('fitD_int', sd=T.std(Y_int), observed=Y_int)

            # Save likelihood
            pm.Deterministic('logp', M.logpt)

        return M

    def sampling(self):
        """This is the sampling that actually runs the model."""
        self.trace = pm.sample(init='advi', model=self.M)

    def fit_ADVI(self):
        with self.M:
            approx = pm.fit(40000, method='fullrank_advi')
            self.trace = approx.sample()

    def profile(self):
        """ Profile the gradient calculation. """
        self.M.profile(pm.theanof.gradient(self.M.logpt, None)).summary()
