import numpy as np, pandas as pds
from .model import getActiveSpecies, runCkineU


class surf_IL2Rb: 
    "Generate values to match the surface IL2Rb measurements used in fitting"
    def __init__(self):
        # times from experiment are hard-coded into this function
        self.ts = np.array([0., 2., 5., 15., 30., 60., 90.])

        # Condense to just IL2Rb
        self.condense = np.zeros(48)
        self.condense[1] = 1

    def calc(self, unkVec):
        "This function uses an unkVec that has the same elements as the unkVec in fit.py"
        
        assert(unkVec.size == 24)
        
        # set IL2 concentrations
        unkVec_IL2_1 = unkVec.copy()
        unkVec_IL2_1[0] = 1.
        
        unkVec_IL2_500 = unkVec.copy()
        unkVec_IL2_500[0] = 500.
        
        # set IL2Ra- values
        unkVecIL2RaMinus_IL2_1 = unkVec_IL2_1.copy()
        unkVecIL2RaMinus_IL2_1[18] = 0.0
        
        unkVecIL2RaMinus_IL2_500 = unkVec_IL2_500.copy()
        unkVecIL2RaMinus_IL2_500[18] = 0.0
        
        # set IL15 concentrations
        unkVec_IL15_1 = unkVec.copy()
        unkVec_IL15_1[1] = 1.
        
        unkVec_IL15_500 = unkVec.copy()
        unkVec_IL15_500[1] = 500.
        
        # set IL2Ra- values
        unkVecIL2RaMinus_IL15_1 = unkVec_IL15_1.copy()
        unkVecIL2RaMinus_IL15_1[18] = 0.0
        
        unkVecIL2RaMinus_IL15_500 = unkVec_IL15_500.copy()
        unkVecIL2RaMinus_IL15_500[18] = 0.0
        
        # calculate IL2 stimulation
        a_yOut, a_retVal = runCkineU(self.ts, unkVec_IL2_1)
        b_yOut, b_retVal = runCkineU(self.ts, unkVec_IL2_500)
        c_yOut, c_retVal = runCkineU(self.ts, unkVecIL2RaMinus_IL2_1)
        d_yOut, d_retVal = runCkineU(self.ts, unkVecIL2RaMinus_IL2_500)
        
        # calculate IL15 stimulation
        e_yOut, e_retVal = runCkineU(self.ts, unkVec_IL15_1)
        f_yOut, f_retVal = runCkineU(self.ts, unkVec_IL15_500)
        g_yOut, g_retVal = runCkineU(self.ts, unkVecIL2RaMinus_IL15_1)
        h_yOut, h_retVal = runCkineU(self.ts, unkVecIL2RaMinus_IL15_500) 
        
        # Condense to just IL2Rb
        condense = np.zeros(48)
        condense[1] = 1
        # species in IL2 family that contain IL2Rb
        condense[4] = 1
        condense[5] = 1
        condense[7] = 1
        condense[8] = 1
        # species in IL15 family that contain IL2Rb
        condense[11] = 1
        condense[12] = 1
        condense[14] = 1
        condense[15] = 1
        
        a = np.dot(a_yOut, condense)
        b = np.dot(b_yOut, condense)
        c = np.dot(c_yOut, condense)
        d = np.dot(d_yOut, condense)
        e = np.dot(e_yOut, condense)
        f = np.dot(f_yOut, condense)
        g = np.dot(g_yOut, condense)
        h = np.dot(h_yOut, condense)

        return (a / a[0], b / b[0], c / c[0], d / d[0], e / e[0], f / f[0], g / g[0], h / h[0])

    
class pstat:
    "Generate values to match the pSTAT5 measurements used in fitting"
    def __init__(self):
        self.cytokC = np.logspace(-3.3, 2.7, 8) # 8 log-spaced values between our two endpoints
    
        npactivity = getActiveSpecies().astype(np.float64)
        self.activity = np.concatenate((npactivity, 0.5*npactivity, np.zeros(4))) # 0.5 is because its the endosome
        self.ts = np.array([500.])
        
    def calc(self, unkVec):
        "This function uses an unkVec that has the same elements as the unkVec in fit.py"
        
        assert(unkVec.size == 24)
        
        # loop over concentrations of IL2
        unkVec_IL2 = np.zeros((24, 8))
        unkVec_IL2_IL2Raminus = unkVec_IL2.copy()
        IL2_yOut = np.ones((8,48))
        IL2_yOut_IL2Raminus = IL2_yOut.copy()
        actVec_IL2 = np.zeros((8))
        actVec_IL2_IL2Raminus = actVec_IL2.copy()
        for ii in range(0,8):
            unkVec_IL2[:, ii] = unkVec.copy()
            unkVec_IL2[0, ii] = self.cytokC[ii]
            
            unkVec_IL2_IL2Raminus[:,ii] = unkVec_IL2[:,ii].copy()
            unkVec_IL2_IL2Raminus[18,ii] = 0.0 # set IL2Ra expression rate to 0
            
            # using Op in hopes of finding time at which activity is maximal and using said time to generate yOut
            IL2_yOut[ii,:], retval_1 = runCkineU(self.ts, unkVec_IL2[:,ii])
            IL2_yOut_IL2Raminus[ii,:], retval_2 = runCkineU(self.ts, unkVec_IL2_IL2Raminus[:,ii])
            
            # dot yOut vectors by activity mask to generate total amount of active species 
            actVec_IL2[ii] = np.dot(IL2_yOut[ii,:], self.activity)
            actVec_IL2_IL2Raminus[ii] = np.dot(IL2_yOut_IL2Raminus[ii,:], self.activity)
            
        print('IL2 stimulation: IL2Ra+')
        print(actVec_IL2)
        print('IL2 stimulation: IL2Ra-')
        print(actVec_IL2_IL2Raminus)
            
        # loop over concentrations of IL15
        unkVec_IL15 = np.zeros((24, 8))
        unkVec_IL15_IL2Raminus = unkVec_IL15.copy()
        IL15_yOut = np.ones((8,48))
        IL15_yOut_IL2Raminus = IL15_yOut.copy()
        actVec_IL15 = np.zeros((8))
        actVec_IL15_IL2Raminus = actVec_IL15.copy()
        for ii in range(0,8):
            unkVec_IL15[:, ii] = unkVec.copy()
            unkVec_IL15[0, ii] = self.cytokC[ii]
            
            unkVec_IL15_IL2Raminus[:,ii] = unkVec_IL15[:,ii].copy()
            unkVec_IL15_IL2Raminus[18,ii] = 0.0 # set IL2Ra expression rate to 0
            
            # using Op in hopes of finding time at which activity is maximal and using said time to generate yOut
            IL15_yOut[ii,:], retval_3 = runCkineU(self.ts, unkVec_IL15[:,ii])
            IL15_yOut_IL2Raminus[ii,:], retval_4 = runCkineU(self.ts, unkVec_IL15_IL2Raminus[:,ii])
            
            # dot yOut vectors by activity mask to generate total amount of active species
            actVec_IL15[ii] = np.dot(IL15_yOut[ii,:], self.activity)
            actVec_IL15_IL2Raminus[ii] = np.dot(IL15_yOut_IL2Raminus[ii,:], self.activity)
            
        print('IL15 stimulation: IL2Ra+')
        print(actVec_IL15)
        print('IL15 stimulation: IL2Ra-')
        print(actVec_IL15_IL2Raminus)

        # Normalize to the maximal activity, put together into one vector
        actVec = np.concatenate((actVec_IL2 / max(actVec_IL2), actVec_IL2_IL2Raminus / max(actVec_IL2_IL2Raminus), actVec_IL15 / max(actVec_IL15), actVec_IL15_IL2Raminus / max(actVec_IL15_IL2Raminus)))
        
        return actVec
    
# def plot(unkVec)