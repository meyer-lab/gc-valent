"""
This file is responsible for performing calculations that allow us to compare our fitting results with the Ring paper in figure1.py
"""
import numpy as np
from .model import getTotalActiveSpecies, runCkineU, runCkineUP, getSurfaceIL2RbSpecies, nParams, getSurfaceGCSpecies


class surf_IL2Rb:
    '''Generate values to match the surface IL2Rb measurements used in fitting'''
    def __init__(self):
        # import function returns from model.py
        self.IL2Rb_species_IDX = getSurfaceIL2RbSpecies()

    def singleCalc(self, unkVec, cytokine, conc, t):
        """ Calculates the surface IL2Rb over time for one condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine, :] = np.ones((unkVec.shape[1])) * conc

        returnn, retVal = runCkineUP(t, np.transpose(unkVec))

        assert retVal >= 0

        a = np.dot(returnn, self.IL2Rb_species_IDX)

        return a

    def calc(self, unkVec, t):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.shape[0] == nParams()

        # set IL2 concentrations
        unkVecIL2RaMinus = unkVec.copy()
        unkVecIL2RaMinus[22, :] = np.zeros((unkVec.shape[1]))
        print("about to call singleCalc for surfIL2Rb")
        print("unkVec.shape: " + str(unkVec.shape))

        # calculate IL2 stimulation
        a = self.singleCalc(unkVec, 0, 1., t)
        b = self.singleCalc(unkVec, 0, 500., t)
        c = self.singleCalc(unkVecIL2RaMinus, 0, 1., t)
        d = self.singleCalc(unkVecIL2RaMinus, 0, 500., t)

        # calculate IL15 stimulation
        e = self.singleCalc(unkVec, 1, 1., t)
        f = self.singleCalc(unkVec, 1, 500., t)
        g = self.singleCalc(unkVecIL2RaMinus, 1, 1., t)
        h = self.singleCalc(unkVecIL2RaMinus, 1, 500., t)

        return np.concatenate((a, b, c, d, e, f, g, h)) / a[0]

class pstat:
    '''Generate values to match the pSTAT5 measurements used in fitting'''
    def __init__(self):
        # import function returns from model.py
        self.activity = getTotalActiveSpecies().astype(np.float64)
        self.ts = np.array([500.]) # was 500. in literature

    def parallelCalc(self, unkVec, cytokine, conc):
        """ Calculates the pSTAT activities in parallel for a 2-D array of unkVec. """
        unkVec = unkVec.copy()
        unkVec[cytokine, :] = np.ones((unkVec.shape[1])) * conc
        print("about to call runCkineUP")

        returnn, retVal = runCkineUP(self.ts, np.transpose(unkVec))
        print("done with runCkineUP")

        assert retVal >= 0

        return np.dot(returnn, self.activity)

    def singleCalc(self, unkVec, cytokine, conc):
        """ Calculates the pSTAT activity for one unkVec condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc

        returnn, retVal = runCkineU(self.ts, unkVec)

        assert retVal >= 0

        return np.dot(returnn, self.activity)
    
    def calc(self, unkVec, cytokC):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''
        assert unkVec.shape[0] == nParams()

        unkVec_IL2Raminus = unkVec.copy()
        unkVec_IL2Raminus[22, :] = np.zeros((unkVec.shape[1])) # set IL2Ra expression rate to 0

        # Calculate activities
        print("about to create actVec_IL2")
        actVec_IL2 = np.fromiter((self.parallelCalc(unkVec, 0, x) for x in cytokC), np.float64)
        actVec_IL2_IL2Raminus = np.fromiter((self.parallelCalc(unkVec_IL2Raminus, 0, x) for x in cytokC), np.float64)
        actVec_IL15 = np.fromiter((self.parallelCalc(unkVec, 1, x) for x in cytokC), np.float64)
        actVec_IL15_IL2Raminus = np.fromiter((self.parallelCalc(unkVec_IL2Raminus, 1, x) for x in cytokC), np.float64)

        # Normalize to the maximal activity, put together into one vector
        actVec = np.concatenate((actVec_IL2, actVec_IL2_IL2Raminus, actVec_IL15, actVec_IL15_IL2Raminus))

        return actVec / np.max(actVec)

class surf_gc:
    """ This class is responsible for calculating the percent of gamma chain on the cell surface. The experimental conditions match those of the surface IL2Rb measurements in Ring et al. """
    def __init__(self):
        # import function returns from model.py
        self.gc_species_IDX = getSurfaceGCSpecies()

    def singleCalc(self, unkVec, cytokine, conc, t):
        """ Calculates the surface gc over time for one condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine, :] = np.ones((unkVec.shape[1])) * conc

        returnn, retVal = runCkineUP(t, np.transpose(unkVec))

        assert retVal >= 0

        a = np.dot(returnn, self.gc_species_IDX)

        return a

    def calc(self, unkVec, t):
        '''This function calls single Calc for all the experimental combinations of interest; it uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.shape[0] == nParams()

        # set IL2 concentrations
        unkVecIL2RaMinus = unkVec.copy()
        unkVecIL2RaMinus[22, :] = np.zeros((unkVec.shape[1]))

        # calculate IL2 stimulation
        a = self.singleCalc(unkVec, 0, 1., t)
        b = self.singleCalc(unkVec, 0, 500., t)
        c = self.singleCalc(unkVecIL2RaMinus, 0, 1., t)
        d = self.singleCalc(unkVecIL2RaMinus, 0, 500., t)

        # calculate IL15 stimulation
        e = self.singleCalc(unkVec, 1, 1., t)
        f = self.singleCalc(unkVec, 1, 500., t)
        g = self.singleCalc(unkVecIL2RaMinus, 1, 1., t)
        h = self.singleCalc(unkVecIL2RaMinus, 1, 500., t)

        return np.concatenate((a, b, c, d, e, f, g, h)) / a[0]
    