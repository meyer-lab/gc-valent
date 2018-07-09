import numpy as np
import matplotlib.pyplot as plt
from .model import getActiveSpecies, runCkineU, getSurfaceIL2RbSpecies


class surf_IL2Rb:
    '''Generate values to match the surface IL2Rb measurements used in fitting'''
    def __init__(self):
        # times from experiment are hard-coded into this function
        self.ts = np.array([0., 2., 5., 15., 30., 60., 90.])

        # Condense to just IL2Rb on surface that is either free, bound to IL2, or bound to IL15
        self.IL2Rb_species_IDX = getSurfaceIL2RbSpecies()

        # percentage value that is used in scaling output
        self.y_max = 10

    def calc(self, unkVec):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.size == 24

        # set IL2 concentrations
        unkVec_IL2_1 = unkVec.copy()
        unkVec_IL2_500 = unkVec.copy()
        unkVec_IL2_1[0], unkVec_IL2_500[0] = 1., 500.

        # set IL2Ra- values
        unkVecIL2RaMinus_IL2_1 = unkVec_IL2_1.copy()
        unkVecIL2RaMinus_IL2_500 = unkVec_IL2_500.copy()
        unkVecIL2RaMinus_IL2_1[18], unkVecIL2RaMinus_IL2_500[18] = 0.0, 0.0

        # set IL15 concentrations
        unkVec_IL15_1 = unkVec.copy()
        unkVec_IL15_500 = unkVec.copy()
        unkVec_IL15_1[1], unkVec_IL15_500[1]  = 1., 500.

        # set IL2Ra- values
        unkVecIL2RaMinus_IL15_1 = unkVec_IL15_1.copy()
        unkVecIL2RaMinus_IL15_500 = unkVec_IL15_500.copy()
        unkVecIL2RaMinus_IL15_1[18], unkVecIL2RaMinus_IL15_500[18] = 0.0, 0.0

        # calculate IL2 stimulation
        a_yOut = runCkineU(self.ts, unkVec_IL2_1)[0]
        b_yOut = runCkineU(self.ts, unkVec_IL2_500)[0]
        c_yOut = runCkineU(self.ts, unkVecIL2RaMinus_IL2_1)[0]
        d_yOut = runCkineU(self.ts, unkVecIL2RaMinus_IL2_500)[0]

        # calculate IL15 stimulation
        e_yOut = runCkineU(self.ts, unkVec_IL15_1)[0]
        f_yOut = runCkineU(self.ts, unkVec_IL15_500)[0]
        g_yOut = runCkineU(self.ts, unkVecIL2RaMinus_IL15_1)[0]
        h_yOut = runCkineU(self.ts, unkVecIL2RaMinus_IL15_500)[0]

        a = np.dot(a_yOut, self.IL2Rb_species_IDX)
        b = np.dot(b_yOut, self.IL2Rb_species_IDX)
        c = np.dot(c_yOut, self.IL2Rb_species_IDX)
        d = np.dot(d_yOut, self.IL2Rb_species_IDX)
        e = np.dot(e_yOut, self.IL2Rb_species_IDX)
        f = np.dot(f_yOut, self.IL2Rb_species_IDX)
        g = np.dot(g_yOut, self.IL2Rb_species_IDX)
        h = np.dot(h_yOut, self.IL2Rb_species_IDX)

        return np.concatenate((a / a[0], b / b[0], c / c[0], d / d[0], e / e[0], f / f[0], g / g[0], h / h[0]))

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

        npactivity = getActiveSpecies().astype(np.float64)
        self.activity = np.concatenate((npactivity, 0.5*npactivity, np.zeros(4))) # 0.5 is because its the endosome
        self.ts = np.array([500.]) # was 500. in literature

        # percentage value that is used in scaling output
        self.y_max = 100

    def calc(self, unkVec):
        '''This function uses an unkVec that has the same elements as the unkVec in fit.py'''

        assert unkVec.size == 24

        # loop over concentrations of IL2
        IL2_yOut = np.ones((self.PTS,48))
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
        IL15_yOut = np.ones((self.PTS,48))
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
