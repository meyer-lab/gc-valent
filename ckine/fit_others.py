"""
This file includes the classes and functions necessary to fit the IL4 and IL7 model to experimental data.
"""
from os.path import join
import theano
import pymc3 as pm, theano.tensor as T, os
import numpy as np, pandas as pds
from .model import getTotalActiveSpecies, nParams, getTotalActiveCytokine
from .differencing_op import runCkineDoseOp, runCkinePreSOp

class IL4_7_activity:
    """ This class is responsible for calculating residuals between model predictions and the data from Gonnord et al. """
    def __init__(self):
        path = os.path.dirname(os.path.abspath(__file__))
        dataIL4 = pds.read_csv(join(path, "./data/Gonnord_S3B.csv")).values # imports IL4 file into pandas array
        dataIL7 = pds.read_csv(join(path, "./data/Gonnord_S3C.csv")).values # imports IL7 file into pandas array

        # units are converted from pg/mL to nM
        self.cytokC_4 = np.array([5., 50., 500., 5000., 50000., 250000.]) / 14900. # 14.9 kDa according to sigma aldrich
        self.cytokC_7 = np.array([1., 10., 100., 1000., 10000., 100000.]) / 17400. # 17.4 kDa according to prospec bio

        self.cytokM = np.zeros((self.cytokC_4.size*2, 6), dtype=np.float64)
        self.cytokM[0:self.cytokC_4.size, 4] = self.cytokC_4
        self.cytokM[self.cytokC_4.size::, 2] = self.cytokC_7

        self.fit_data = np.concatenate((dataIL4[:, 1], dataIL4[:, 2], dataIL7[:, 1], dataIL7[:, 2])) # measurements aren't normalized

    def calc(self, unkVec, scales):
        """ Simulate the experiment with different ligand stimulations and compare with experimental data. """
        Op = runCkineDoseOp(tt=np.array(10.), condense=getTotalActiveSpecies().astype(np.float64), conditions=self.cytokM)

        # Run the experiment
        outt = Op(unkVec)

        actVecIL4 = outt[0:self.cytokC_4.size]
        actVecIL7 = outt[self.cytokC_4.size:self.cytokC_4.size*2]

        # Multiply by scaling constants and put together in one vector
        actVec = T.concatenate((actVecIL4 * scales[0], actVecIL4 * scales[0], actVecIL7 * scales[1], actVecIL7 * scales[1]))

        # return residual
        return self.fit_data - actVec

class crosstalk:
    """ This class performs the calculations necessary in order to fit our model to Gonnord Fig S3D. """
    def __init__(self):
        self.activity = getTotalActiveSpecies().astype(np.float64)
        self.ts = np.array([10.]) # was 10. in literature
        self.IL4_stim_conc = 100. / 14900. # concentration used for IL4 stimulation
        self.IL7_stim_conc = 50. / 17400. # concentration used for IL7 stimulation
        self.cytokM = np.zeros((2, 6), dtype=np.float64)
        self.cytokM[0, 4] = self.IL4_stim_conc
        self.cytokM[1, 2] = self.IL7_stim_conc
        
        path = os.path.dirname(os.path.abspath(__file__))
        data = pds.read_csv(join(path, "./data/Gonnord_S3D.csv")).values
        self.fit_data = np.concatenate((data[:, 1], data[:, 2], data[:, 3], data[:, 6], data[:, 7], data[:, 8]))
        self.pre_IL7 = data[:, 0]   # concentrations of IL7 used as pretreatment
        self.pre_IL4 = data[:, 5]   # concentrations of IL4 used as pretreatment


    def singleCalc(self, unkVec, pre_cytokine, pre_conc, stim_cytokine, stim_conc):
        """ This function generates the active vector for a given unkVec, cytokine used for inhibition and concentration of pretreatment cytokine. """
        unkVec2 = T.set_subtensor(unkVec[pre_cytokine], pre_conc)
        ligands = np.zeros((6))
        ligands[stim_cytokine] = stim_conc
        ligands[pre_cytokine] = pre_conc

        Op = runCkinePreSOp(tpre=self.ts, ts=self.ts, postlig=ligands)

        # perform the experiment
        outt = Op(unkVec2)

        return getTotalActiveCytokine(stim_cytokine, outt) # only look at active species associated with stimulation cytokine

    def singleCalc_no_pre(self, unkVec):
        ''' This function generates the active vector for a given unkVec, cytokine, and concentration. '''
        Op = runCkineDoseOp(tt=self.ts, condense=getTotalActiveSpecies().astype(np.float64), conditions=self.cytokM)

        # Run the experiment
        outt = Op(unkVec)

        actVecIL4 = outt[0]
        actVecIL7 = outt[1]
        return actVecIL4, actVecIL7


    def calc(self, unkVec):
        """ Generates residual calculation that compares model to data. """
        # with no pretreatment
        IL4stim_no_pre, IL7stim_no_pre = self.singleCalc_no_pre(unkVec)

        # add 6 zeros to front of unkVec to handle ligands
        unkVec = T.concatenate((T.zeros(6, dtype=np.float64), unkVec))

        # IL7 pretreatment with IL4 stimulation
        actVec_IL4stim = T.stack((list(self.singleCalc(unkVec, 2, x, 4, self.IL4_stim_conc) for x in self.pre_IL7)))

        # IL4 pretreatment with IL7 stimulation
        actVec_IL7stim = T.stack((list(self.singleCalc(unkVec, 4, x, 2, self.IL7_stim_conc) for x in self.pre_IL4)))

        case1 = 1-(actVec_IL4stim/IL4stim_no_pre) * 100.    # % inhibition of IL4 act. after IL7 pre.
        case2 = 1-(actVec_IL7stim/IL7stim_no_pre) * 100.    # % inhibition of IL7 act. after IL4 pre.
        inh_vec = T.concatenate((case1, case1, case1, case2, case2, case2 ))   # mimic order of CSV file

        return inh_vec - self.fit_data


class build_model:
    """Going to load the data from the CSV file at the very beginning of when build_model is called... needs to be separate member function to avoid uploading file thousands of times."""
    def __init__(self):
        self.act = IL4_7_activity()
        self.cross = crosstalk()
        self.M = self.build()

    def build(self, pretreat=False):
        """The PyMC model that incorporates Bayesian Statistics in order to store what the likelihood of the model is for a given point."""
        M = pm.Model()

        with M:
            kfwd = pm.Lognormal('kfwd', mu=np.log(0.00001), sd=1, shape=1)
            nullRates = T.ones(6, dtype=np.float64) # associated with IL2 and IL15
            Tone = T.ones(1, dtype=np.float64)
            Tzero = T.zeros(1, dtype=np.float64)
            k27rev = pm.Lognormal('k27rev', mu=np.log(0.1), sd=10, shape=1) # associated with IL7
            k33rev = pm.Lognormal('k33rev', mu=np.log(0.1), sd=10, shape=1) # associated with IL4
            endo_activeEndo = pm.Lognormal('endo', mu=np.log(0.1), sd=0.1, shape=2)
            sortF = pm.Beta('sortF', alpha=20, beta=40, testval=0.333, shape=1)*0.95
            kRec_kDeg = pm.Lognormal('kRec_kDeg', mu=np.log(0.1), sd=0.1, shape=2)
            GCexpr = (328. * endo_activeEndo[0]) / (1. + ((kRec_kDeg[0]*(1.-sortF)) / (kRec_kDeg[1]*sortF))) # constant according to measured number per cell
            IL7Raexpr = (2591. * endo_activeEndo[0]) / (1. + ((kRec_kDeg[0]*(1.-sortF)) / (kRec_kDeg[1]*sortF))) # constant according to measured number per cell
            IL4Raexpr = (254. * endo_activeEndo[0]) / (1. + ((kRec_kDeg[0]*(1.-sortF)) / (kRec_kDeg[1]*sortF))) # constant according to measured number per cell
            scales = pm.Lognormal('scales', mu=np.log(100), sd=np.log(25), shape=2) # create scaling constants for activity measurements

            unkVec = T.concatenate((kfwd, nullRates, k27rev, Tone, k33rev, Tone, endo_activeEndo, sortF, kRec_kDeg))
            unkVec = T.concatenate((unkVec, Tzero, Tzero, GCexpr, Tzero, IL7Raexpr, Tzero, IL4Raexpr, Tzero)) # indexing same as in model.hpp

            Y_int = self.act.calc(unkVec, scales) # fitting the data based on act.calc for the given parameters

            pm.Deterministic('Y_int', T.sum(T.square(Y_int)))

            pm.Normal('fitD_int', sd=700, observed=Y_int)
            
            if pretreat == True:
                print('inside if-statement')
                Y_cross = self.cross.calc(unkVec)   # fitting the data based on cross.calc
                pm.Deterministic('Y_cross', T.sum(T.square(Y_cross)))
                pm.Normal('fitD_cross', sd=0.1, observed=Y_cross)

            # Save likelihood
            pm.Deterministic('logp', M.logpt)

        return M

    def sampling(self):
        """This is the sampling that actually runs the model."""
        self.trace = pm.sample(init='advi', model=self.M)

    def fit_ADVI(self):
        """ Running fit_advi instead of true sampling. """
        with self.M:
            approx = pm.fit(40000, method='fullrank_advi')
            self.trace = approx.sample()

    def profile(self):
        """ Profile the gradient calculation. """
        self.M.profile(pm.theanof.gradient(self.M.logpt, None)).summary()

