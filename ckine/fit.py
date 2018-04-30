"""
This file includes the classes and functions necessary to fit the IL2 model to the experimental data.
"""
import pymc3 as pm, theano.tensor as T, os
from theano import shared
import numpy as np, pandas as pds
from .model import runCkineU, getActiveSpecies
from .differencing_op import centralDiff, runCkineOp


def surf_IL2Rb(rxntraf, IL2_conc):
    # times from experiment are hard-coded into this function
    ts = np.array(([0.01, 2., 5., 15., 30., 60., 90.]))

    rxntraf[0] = IL2_conc # the concentration of IL2 is rxnRates[0]

    ys, retVal = runCkineU(ts, rxntraf)

    if retVal < 0:
        return -100

    return 10. * (ys[:, 1] / ys[0, 1]) # % sIL2Rb relative to initial amount


class IL2Rb_trafficking:
    def __init__(self):
        path = os.path.dirname(os.path.abspath(__file__))
        data = pds.read_csv(os.path.join(path, 'data/IL2Ra+_surface_IL2RB_datasets.csv')) # imports csv file into pandas array
        self.numpy_data = data.as_matrix() # all of the IL2Rb trafficking data with IL2Ra+... first row contains headers... 9 columns and 8 rows... first column is time
        data2 = pds.read_csv(os.path.join(path, "data/IL2Ra-_surface_IL2RB_datasets.csv"))
        self.numpy_data2 = data2.as_matrix() # all of the IL2Rb trafficking data with IL2Ra-... first row contains headers... 9 columns and 8 rows... first column is time

        self.concs = 14

    def calc(self, tfR):
        # IL2Ra- cells
        tfR2 = tfR.copy()
        tfR2[15] = 0.0

        diff1 = surf_IL2Rb(tfR, 1) - self.numpy_data[:, 1] # the second column of numpy_data has all the 1nM IL2Ra+ data
        diff2 = surf_IL2Rb(tfR, 500) - self.numpy_data[:, 5] # the sixth column of numpy_data has all the 500 nM IL2Ra+ data
        diff3 = surf_IL2Rb(tfR2, 1) - self.numpy_data2[:, 1] # the second column of numpy_data2 has all the 1nM IL2Ra- data
        diff4 = surf_IL2Rb(tfR2, 500) - self.numpy_data2[:, 5] # the sixth column of numpy_data2 has all the 500 nM IL2Ra- data

        all_diffs = np.concatenate((diff1, diff2, diff3, diff4))

        return all_diffs
    
# this takes all the desired IL2 values we want to test and gives us the maximum activity value
# IL2 values pretty much ranged from 5 x 10**-4 to 500 nm with 8 points in between
# need the theano decorator to get around the fact that there are if-else statements when running odeint but
# we don't necessarily know the values for the rxn rates when we call our model
class IL2_15_activity:
    def __init__(self):
        """This loads the experiment data and saves it as a member matrix and it also makes a vector of the IL15 concentrations that we are going to take care of."""
        path = os.path.dirname(os.path.abspath(__file__))
        data = pds.read_csv(os.path.join(path, "./data/IL2_IL15_extracted_data.csv")).as_matrix() # imports csv file into pandas array
        dataIL2 = pds.read_csv(os.path.join(path, "./data/IL2_IL15_extracted_data.csv")).as_matrix() # imports csv file into pandas array
        self.cytokC = np.logspace(-3.3, 2.7, 8) # 8 log-spaced values between our two endpoints
        self.fit_data = np.concatenate((data[:, 7], data[:, 3], dataIL2[:, 6], dataIL2[:, 2])) #the IL15_IL2Ra- data is within the 4th column (index 3)
        # the IL2_IL2Ra- data is within the 3rd column (index 2)
        npactivity = getActiveSpecies().astype(np.float64)
        self.activity = shared(np.concatenate((npactivity, 0.5*npactivity, np.zeros(4)))) # 0.5 is because its the endosome

    def calc(self, unkVec):
        """Simulate the experiment with IL15. It is making a list of promises which will be calculated and returned as output."""
        # Convert the vector of values to dicts

        # IL2Ra- cells have same IL15 activity, so we can just reuse same solution
        Op = runCkineOp(ts=np.array(500.))

        # Loop over concentrations of IL15
        actVec = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[1], x))), self.cytokC))) # Change condensation here for activity

        # Loop over concentrations of IL2
        actVecIL2 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[0], x))), self.cytokC)))

        unkVecIL2RaMinus = T.set_subtensor(unkVec[15], 0.0) # Set IL2Ra to zero

        # Loop over concentrations of IL2, IL2Ra-/-
        actVecIL2RaMinus = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVecIL2RaMinus[0], x))), self.cytokC)))

        # Normalize to the maximal activity, put together into one vector
        actVec = T.concatenate((actVec / T.max(actVec), actVec / T.max(actVec), actVecIL2 / T.max(actVecIL2), actVecIL2RaMinus / T.max(actVecIL2RaMinus)))
        # value we're trying to minimize is the distance between the y-values on points of the graph that correspond to the same IL2 values
        return self.fit_data - actVec
    

class build_model:
    """Going to load the data from the CSV file at the very beginning of when build_model is called... needs to be separate member function to avoid uploading file thousands of times."""
    def __init__(self):
        self.dst15 = IL2_15_activity()
        self.IL2Rb = IL2Rb_trafficking()
        self.M = self.build()

    def build(self):
        """The PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point."""
        M = pm.Model()

        with M:
            rxnrates = pm.Lognormal('rxn', mu=np.log(0.1), sd=1., shape=8) # first 3 are IL2, second 5 are IL15, kfwd is first element (used in both 2&15)
            endo_activeEndo = pm.Lognormal('endo', mu=np.log(0.1), sd=1., shape=2)
            kRec_kDeg = pm.Lognormal('kRec_kDeg', mu=np.log(0.1), sd=1., shape=2)
            Rexpr = pm.Lognormal('IL2Raexpr', sd=1., shape=4) # Expression: IL2Ra, IL2Rb, gc, IL15Ra
            sortF = pm.Beta('sortF', alpha=2, beta=7, testval=0.1)

            ligands = T.zeros(4, dtype=np.float64)

            unkVec = T.concatenate((ligands, rxnrates, T.zeros(3, dtype=np.float64),
                                    endo_activeEndo, T.stack(sortF), kRec_kDeg, Rexpr, T.zeros(2, dtype=np.float64)))

            Y_15 = self.dst15.calc(unkVec) # fitting the data based on dst15.calc for the given parameters
            Y_int = centralDiff(self.IL2Rb)(unkVec) # fitting the data based on dst.calc for the given parameters

            pm.Deterministic('Y_15', T.sum(T.square(Y_15)))
            pm.Deterministic('Y_int', T.sum(T.square(Y_int)))

            pm.Normal('fitD_15', sd=T.std(Y_15), observed=Y_15)
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
