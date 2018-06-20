import theano.tensor as T, os
from os.path import join
from theano import shared
import numpy as np, pandas as pds
from .model import getActiveSpecies, runCkineU
from .differencing_op import runCkineOp, runCkineKineticOp

class surf_IL2Rb: 
    "the only difference between calc_1 and calc_2 is the values in the unkVec that are passed in"
    def __init__(self):
        # times from experiment are hard-coded into this function
        self.ts = np.array([0., 2., 5., 15., 30., 60., 90.])

        # Condense to just IL2Rb
        self.condense = np.zeros(48)
        self.condense[1] = 1

    def calc_1(self):
        "This function uses an unkVec from fitting surface IL2Rb data only. This data should closely match the experimental figures"
        
        # posteriors taken from line 170 of fitting_IL2Rb_only.csv
        # NOTE: these results were fit to IL2 stimulation only
        unkVec = np.array([0., 0., 0., 0., 8.33E-05, 0.115818087, 0.128608864, 0.089890767, 0.098851823, 0.106913961, 0.107722847, 0.117471386, 0.114016182, 0.095431878, 0.101391307, 0.340647077, 0.090746998, 0.097587796, 1.013323342, 0.925694209, 0.918209005, 0.958374837, 0., 0.])

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

    def calc_2(self):
        "This function uses an unkVec from fitting surface IL2Rb and IL2 pSTAT data. This fitting process had high error so the return values mightly highly differ from the experimental figures "
        
        # posteriors taken from line 270 of full-fitting-no-IL15-hard-bound.csv
        unkVec = np.array([0., 0., 0., 0., 0.00012865, 0.158037088, 0.285547466, 0.102937877, 0.399706229, 0.111318564, 0.061811833, 0.042950052, 0.283564271, 0.088275756, 0.08486402, 0.399519643, 0.05414861, 0.253249545, 0.509200117, 0.537992752, 0.667527049, 3.852037782, 0., 0.])
        
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
    def __init__(self):
        self.cytokC = np.logspace(-3.3, 2.7, 8) # 8 log-spaced values between our two endpoints
    
        npactivity = getActiveSpecies().astype(np.float64)
        self.activity = np.concatenate((npactivity, 0.5*npactivity, np.zeros(4))) # 0.5 is because its the endosome
        self.ts = np.array(500.)
        
    def calc_1(self):
        # posteriors taken from line 270 of full-fitting-no-IL15-hard-bound.csv
        # order of elements reagganged to match unkVec in fit.py
        vec = np.array([0., 0., 0., 0., 0.00012865, 0.158037088, 0.285547466, 0.102937877, 0.399706229, 0.111318564, 0.061811833, 0.042950052, 0.283564271, 0.088275756, 0.08486402, 0.399519643, 0.05414861, 0.253249545, 0.509200117, 0.537992752, 0.667527049, 3.852037782, 0., 0.])
        
        Op = runCkineOp(ts=np.array(500.))
        
        # loop over concentrations of IL2
        unkVec_IL2 = np.zeros((24, 8))
        unkVec_IL2_IL2Raminus = unkVec_IL2.copy()
        IL2_yOut = np.ones((8,48))
        IL2_yOut_IL2Raminus = IL2_yOut.copy()
        act_IL2 = np.ones((8))
        act_IL2_IL2Raminus = act_IL2.copy()
        for ii in range(0,8):
            unkVec_IL2[:, ii] = vec.copy()
            unkVec_IL2[0, ii] = self.cytokC[ii]
            
            unkVec_IL2_IL2Raminus[:,ii] = unkVec_IL2[:,ii].copy()
            unkVec_IL2_IL2Raminus[18,ii] = 0.0 # set IL2Ra expression rate to 0
            
            # using Op in hopes of finding time at which activity is maximal and using said time to generate yOut
            IL2_yOut[ii,:] = Op(unkVec_IL2[:,ii])
            IL2_yOut_IL2Raminus[ii,:] = Op(unkVec_IL2_Il2Raminus[:,ii])
            
            # dot yOut vectors by activity mask to generate total amount of active species
            act_IL2[ii] = np.dot(IL2_yOut[ii,:], self.activity)
            act_IL2_IL2Raminus[ii] = np.dot(IL2_yOut_IL2Raminus[ii,:], self.activity)
            
        print(act_IL2)
        print(act_IL2_IL2Raminus)
            
        # loop over concentrations of IL15
        unkVec_IL15 = np.zeros((24, 8))
        unkVec_IL15_IL2Raminus = unkVec_IL15.copy()
        IL15_yOut = np.ones((8,48))
        IL15_yOut_IL2Raminus = IL15_yOut.copy()
        act_IL15 = np.ones((8))
        act_IL15_IL2Raminus = act_IL15.copy()
        for ii in range(0,8):
            unkVec_IL15[:, ii] = vec.copy()
            unkVec_IL15[1, ii] = self.cytokC[ii]
            
            unkVec_IL15_IL2Raminus[:,ii] = unkVec_IL15[:,ii].copy()
            unkVec_IL15_IL2Raminus[18,ii] = 0.0 # set IL2Ra expression rate to 0
            
            # using Op in hopes of finding time at which activity is maximal and using said time to generate yOut
            IL15_yOut[ii,:] = Op(unkVec_IL15[:,ii])
            IL15_yOut_IL2Raminus[ii,:] = Op(unkVec_IL15_Il2Raminus[:,ii])
            
            # dot yOut vectors by activity mask to generate total amount of active species
            act_IL15[ii] = np.dot(IL15_yOut[ii,:], self.activity)
            act_IL15_IL2Raminus[ii] = np.dot(IL15_yOut_IL2Raminus[ii,:], self.activity)
            
        print(act_IL15)
        print(act_IL15_IL2Raminus)

        # Normalize to the maximal activity, put together into one vector
        # actVec = T.concatenate((actVec / T.max(actVec), actVec / T.max(actVec), actVecIL2 / T.max(actVecIL2), actVecIL2RaMinus / T.max(actVecIL2RaMinus)))
        actVec = np.concatenate((act_IL2 / max(act_IL2), act_IL2_IL2Raminus / max(act_IL2_IL2Raminus), act_IL15 / max(act_IL15), act_IL15_IL2Raminus / max(act_IL15_IL2Raminus)))
        # value we're trying to minimize is the distance between the y-values on points of the graph that correspond to the same IL2 values
        return actVec