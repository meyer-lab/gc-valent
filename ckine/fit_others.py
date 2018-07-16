"""
This file includes the classes and functions necessary to fit the IL4, IL7 and IL21 model to the experimental data.
"""
import pymc3 as pm, theano.tensor as T, os
from os.path import join
from theano import shared
import numpy as np, pandas as pds
from .model import getTotalActiveSpecies
from .differencing_op import runCkineOp, runCkineKineticOp

class IL4_7_activity:
    def __init__(self):
        """This loads the experiment data and saves it as a member matrix and it also makes a vector of the ligand concentrations that we are going to take care of."""
        path = os.path.dirname(os.path.abspath(__file__))
        dataIL4 = pds.read_csv(join(path, "./data/Gonnard_S3B.csv")).as_matrix() # imports IL4 file into pandas array
        dataIL7 = pds.read_csv(join(path, "./data/Gonnard_S3C.csv")).as_matrix() # imports IL21 file into pandas array
        
        # units are converted from pg/mL to nM
        self.cytokC_4 = np.array([5., 50., 500., 5000., 50000., 250000.]) / 14900. # 14.9 kDa according to sigma aldrich
        self.cytokC_7 = np.array([1., 10., 100., 1000., 10000., 100000.]) / 17400. # 17.4 kDa according to prospec bio

        self.fit_data = np.concatenate((dataIL4[:, 1], dataIL4[:, 2], dataIL7[:, 1], dataIL7[:, 2])) # the measurements are not normalized
        self.activity = getTotalActiveSpecies().astype(np.float64)


    def calc(self, unkVec):
        """Simulate the experiment with different ligand stimulations. It is making a list of promises which will be calculated and returned as output."""
        Op = runCkineOp(ts=np.array(10.))

        # Loop over concentrations of IL4
        actVecIL4 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[4], x))), self.cytokC_4)))

        # Loop over concentrations of IL7
        actVecIL7 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[2], x))), self.cytokC_7)))

        # Normalize to the scaling constants, put together into one vector
        # TODO: make sure indexing in unkVec is correct for scales
        actVec = T.concatenate((actVecIL4 * unkVec[-2], actVecIL4 * unkVec[-2], actVecIL7 * unkVec[-1], actVecIL7 * unkVec[-1]))

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
            GCexpr = (328. * endo_activeEndo[0]) / (1 + (kRec_kDeg[0] / kRec_kDeg[1])) # constant according to measured number per cell
            IL7Raexpr = (2591. * endo_activeEndo[0]) / (1 + (kRec_kDeg[0] / kRec_kDeg[1])) # constant according to measured number per cell
            IL4Raexpr = (254. * endo_activeEndo[0]) / (1 + (kRec_kDeg[0] / kRec_kDeg[1])) # constant according to measured number per cell
            IL21Raexpr = pm.Lognormal('IL21Raexpr', sd=0.1, shape=1)
            sortF = pm.Beta('sortF', alpha=20, beta=40, testval=0.333)*0.95
            # TODO: double check the priors for scales seem reasonable
            scales = pm.Lognormal('scales', mu=np.log(10), sd=np.log(25), shape=2) # create scaling constants for activity measurements

            ligands = T.zeros(6, dtype=np.float64)
            
            # TODO: make sure three measured values are inputted correctly

            unkVec = T.concatenate((ligands, T.stack(kfwd), rxnrates, endo_activeEndo, T.stack(sortF), kRec_kDeg, T.zeros(2, dtype=np.float64), GCexpr, T.zeros(1, dtype=np.float64), IL7Raexpr, T.zeros(1, dtype=np.float64), IL4Raexpr, T.stack(IL21Raexpr), scales)) # receptor expression indexing same as in model.cpp

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
