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

        # Loop over concentrations of IL15
        actVec = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[1], x))), self.cytokC))) # Change condensation here for activity

        # Loop over concentrations of IL2
        actVecIL2 = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVec[0], x))), self.cytokC)))

        unkVecIL2RaMinus = T.set_subtensor(unkVec[18], 0.0) # Set IL2Ra to zero

        # Loop over concentrations of IL2, IL2Ra-/-
        actVecIL2RaMinus = T.stack(list(map(lambda x: T.dot(self.activity, Op(T.set_subtensor(unkVecIL2RaMinus[0], x))), self.cytokC)))

        # Normalize to the maximal activity, put together into one vector
        actVec = T.concatenate((actVec / T.max(actVec), actVec / T.max(actVec), actVecIL2 / T.max(actVecIL2), actVecIL2RaMinus / T.max(actVecIL2RaMinus)))

        # value we're trying to minimize is the distance between the y-values on points of the graph that correspond to the same IL2 values
        return self.fit_data - actVec
