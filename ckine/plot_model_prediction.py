import numpy as np
import matplotlib.pyplot as plt
from .model import getTotalActiveSpecies, runCkineU, getSurfaceIL2RbSpecies, nSpecies, nParams


class surf_IL2Rb:
    '''Generate values to match the surface IL2Rb measurements used in fitting'''
    def __init__(self):
        # times from experiment are hard-coded into this function
        self.ts = np.array([0., 2., 5., 15., 30., 60., 90.])

        # import function returns from model.py
        self.IL2Rb_species_IDX = getSurfaceIL2RbSpecies()
        self.nParams = nParams()

        # percentage value that is used in scaling output
        self.y_max = 10

    def singleCalc(self, unkVec, cytokine, conc):
        """ Calculates the surface IL2Rb over time for one condition. """
        unkVec = unkVec.copy()
        unkVec[cytokine] = conc

        returnn, retVal = runCkineU(self.ts, unkVec)

        assert retVal >= 0

        a = np.dot(returnn, self.IL2Rb_species_IDX)

        return a / a[0]

    def calc(self, unkVec):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.size == self.nParams

        # set IL2 concentrations
        unkVecIL2RaMinus = unkVec.copy()
        unkVecIL2RaMinus[18] = 0.

        # calculate IL2 stimulation
        a = self.singleCalc(unkVec, 0, 1.)
        b = self.singleCalc(unkVec, 0, 500.)
        c = self.singleCalc(unkVecIL2RaMinus, 0, 1.)
        d = self.singleCalc(unkVecIL2RaMinus, 0, 500.)

        # calculate IL15 stimulation
        e = self.singleCalc(unkVec, 1, 1.)
        f = self.singleCalc(unkVec, 1, 500.)
        g = self.singleCalc(unkVecIL2RaMinus, 1, 1.)
        h = self.singleCalc(unkVecIL2RaMinus, 1, 500.)

        return np.concatenate((a, b, c, d, e, f, g, h))

    def plot_structure(self, IL2vec, IL15vec, title):
        plt.title(title)
        plt.scatter(self.ts, IL2vec, color='r', label='IL2', alpha=0.7)
        plt.scatter(self.ts, IL15vec, color='g', label='IL15', alpha=0.7)
        # plt.ylim(0,(y_max + (0.2 * y_max)))
        plt.ylabel("Surface IL2Rb (% x " + str(self.y_max) + ')')
        plt.xlabel("Time (min)")
        plt.legend()
        plt.show()

    def plot(self, unkVec):
        output = self.calc(unkVec) * self.y_max
        IL2_1_plus = output[0:7]
        IL2_500_plus = output[7:14]
        IL2_1_minus = output[14:21]
        IL2_500_minus = output[21:28]
        IL15_1_plus = output[28:35]
        IL15_500_plus = output[35:42]
        IL15_1_minus = output[42:49]
        IL15_500_minus = output[49:56]

        self.plot_structure(IL2_1_minus, IL15_1_minus, '1 nM and IL2Ra-')
        self.plot_structure(IL2_500_minus, IL15_500_minus, "500 nM and IL2Ra-")
        self.plot_structure(IL2_1_plus, IL15_1_plus, "1 nM and IL2Ra+")
        self.plot_structure(IL2_500_plus, IL15_500_plus, "500 nM and IL2Ra+")


class pstat:
    '''Generate values to match the pSTAT5 measurements used in fitting'''
    def __init__(self):
        self.PTS = 25
        self.cytokC = np.logspace(-3.3, 2.7, self.PTS) # 8 log-spaced values between our two endpoints
        
        # import function returns from model.py
        self.nParams = nParams()
        self.nSpecies = nSpecies()
        self.activity = getTotalActiveSpecies().astype(np.float64)
        
        self.ts = np.array([500.]) # was 500. in literature

        # percentage value that is used in scaling output
        self.y_max = 100

    def calc(self, unkVec):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.size == self.nParams

        # loop over concentrations of IL2
        IL2_yOut = np.ones((self.PTS, self.nSpecies))
        IL2_yOut_IL2Raminus = IL2_yOut.copy()
        actVec_IL2 = np.zeros((self.PTS))
        actVec_IL2_IL2Raminus = actVec_IL2.copy()
        for ii in range(0,self.PTS):
            unkVec_IL2 = unkVec.copy()
            unkVec_IL2[0] = self.cytokC[ii]

            unkVec_IL2_IL2Raminus = unkVec_IL2.copy()
            unkVec_IL2_IL2Raminus[18] = 0.0 # set IL2Ra expression rate to 0

            IL2_yOut[ii,:], retval_1 = runCkineU(self.ts, unkVec_IL2)
            IL2_yOut_IL2Raminus[ii,:], retval_2 = runCkineU(self.ts, unkVec_IL2_IL2Raminus)

            if retval_1 < 0:
                print("runCkineU failed for IL2 stimulation in IL2Ra+ cells")
                print("failure occured during iteration " + str(ii))

            if retval_2 < 0:
                print("runCkineU failed for IL2 stimulation in IL2Ra- cells")
                print("failure occured during iteration " + str(ii))

            # dot yOut vectors by activity mask to generate total amount of active species
            actVec_IL2[ii] = np.dot(IL2_yOut[ii,:], self.activity)
            actVec_IL2_IL2Raminus[ii] = np.dot(IL2_yOut_IL2Raminus[ii,:], self.activity)

        # loop over concentrations of IL15
        IL15_yOut = np.ones((self.PTS, self.nSpecies))
        IL15_yOut_IL2Raminus = IL15_yOut.copy()
        actVec_IL15 = np.zeros((self.PTS))
        actVec_IL15_IL2Raminus = actVec_IL15.copy()
        for ii in range(0,self.PTS):
            unkVec_IL15 = unkVec.copy()
            unkVec_IL15[1] = self.cytokC[ii]

            unkVec_IL15_IL2Raminus = unkVec_IL15.copy()
            unkVec_IL15_IL2Raminus[18] = 0.0 # set IL2Ra expression rate to 0

            IL15_yOut[ii,:], retval_3 = runCkineU(self.ts, unkVec_IL15)
            IL15_yOut_IL2Raminus[ii,:], retval_4 = runCkineU(self.ts, unkVec_IL15_IL2Raminus)

            if retval_3 < 0:
                print("runCkineU failed for IL15 stimulation in IL2Ra+ cells")
                print("failure occured during iteration " + str(ii))

            if retval_4 < 0:
                print("runCkineU failed for IL15 stimulation in IL2Ra- cells")
                print("failure occured during iteration " + str(ii))

            # dot yOut vectors by activity mask to generate total amount of active species
            actVec_IL15[ii] = np.dot(IL15_yOut[ii,:], self.activity)
            actVec_IL15_IL2Raminus[ii] = np.dot(IL15_yOut_IL2Raminus[ii,:], self.activity)

        # Normalize to the maximal activity, put together into one vector
        actVec = np.concatenate((actVec_IL2 / np.max(actVec_IL2), actVec_IL2_IL2Raminus / np.max(actVec_IL2_IL2Raminus), actVec_IL15 / np.max(actVec_IL15), actVec_IL15_IL2Raminus / np.max(actVec_IL15_IL2Raminus)))

        return actVec

    def plot_structure(self, IL2vec, IL15vec, title):
        plt.title(title)
        plt.scatter(np.log10(self.cytokC), IL2vec, color='r', alpha=0.5, label="IL2")
        plt.scatter(np.log10(self.cytokC), IL15vec, color='g', alpha=0.5, label='IL15')
        plt.ylim(0,(self.y_max + (0.25*self.y_max)))
        plt.ylabel('Maximal p-STAT5 (% x ' + str(self.y_max) + ')')
        plt.xlabel('log10 of cytokine concentration (nM)')
        plt.legend()
        plt.show()

    def plot(self, unkVec):
        output = self.calc(unkVec) * self.y_max
        IL2_plus = output[0:self.PTS]
        IL2_minus = output[self.PTS:(self.PTS*2)]
        IL15_plus = output[(self.PTS*2):(self.PTS*3)]
        IL15_minus = output[(self.PTS*3):(self.PTS*4)]

        self.plot_structure(IL2_minus, IL15_minus, "IL2Ra- YT-1 cells")
        self.plot_structure(IL2_plus, IL15_plus, "IL2Ra+ YT-1 cells")
