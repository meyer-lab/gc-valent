import pymc3 as pm
import theano.tensor as T
import numpy as np
from .model import getTotalActiveSpecies, receptor_expression
from .differencing_op import runCkineDoseOp
from .imports import import_Rexpr, import_pstat


def commonTraf():
    """ Set the common trafficking parameter priors. """
    kfwd = pm.Lognormal('kfwd', mu=np.log(0.001), sd=0.5, shape=1)
    endo = pm.Lognormal('endo', mu=np.log(0.1), sd=0.1, shape=1)
    activeEndo = pm.Lognormal('activeEndo', sd=0.1, shape=1)
    kRec = pm.Lognormal('kRec', mu=np.log(0.1), sd=0.1, shape=1)
    kDeg = pm.Lognormal('kDeg', mu=np.log(0.01), sd=0.2, shape=1)
    sortF = pm.Beta('sortF', alpha=12, beta=80, shape=1)
    return kfwd, endo, activeEndo, kRec, kDeg, sortF


class pSTAT_activity:
    """ Calculating the pSTAT activity residuals for IL2 and IL15 stimulation in distinct cell lines from Visterra. """
    
    def __init__(self):
        self.ts = np.array([0.5, 1., 2., 4.]) * 60.
        # might have to move these lines into model class and pass the relevant data as function args
        _, self.receptor_data, self.cell_names_receptor = import_Rexpr()
        self.cytokC, self.cell_names_pstat, self.IL2_data, self.IL15_data = import_pstat()

        self.cytokM = np.zeros((self.cytokC.size * 2, 6), dtype=np.float64)
        self.cytokM[0:self.cytokC.size, 0] = self.cytokC
        self.cytokM[self.cytokC.size::, 1] = self.cytokC
        
    def calc(self, unkVec, scale):
        """ Simulate the STAT5 measurements and return residuals between model prediction and experimental data. """

        # Find amount of active species for each simulation
        Op = runCkineDoseOp(tt=self.ts, condense=getTotalActiveSpecies().astype(np.float64), conditions=self.cytokM)
        
        tot_res = T.zeros((self.receptor_data.shape[0], self.ts.size*2, self.cytokC.size))
        for i, name in enumerate(self.cell_names_pstat):
            if name == "CD8+" or name == "Naive Treg":
                print("caught", name)
                # plot matching experimental and predictive pSTAT data for the same cell type
                for j in range(self.receptor_data.shape[0]):
                    if self.cell_names_pstat[i] == self.cell_names_receptor[j]:
                        unkVec = T.set_subtensor(unkVec[16], receptor_expression(self.receptor_data[j, 0], unkVec[11], unkVec[14], unkVec[13], unkVec[15]))  # IL2Ra
                        unkVec = T.set_subtensor(unkVec[17], receptor_expression(self.receptor_data[j, 1], unkVec[11], unkVec[14], unkVec[13], unkVec[15]))  # IL2Rb
                        unkVec = T.set_subtensor(unkVec[18], receptor_expression(self.receptor_data[j, 2], unkVec[11], unkVec[14], unkVec[13], unkVec[15]))  # gc
                        unkVec = T.set_subtensor(unkVec[19], receptor_expression(self.receptor_data[j, 3], unkVec[11], unkVec[14], unkVec[13], unkVec[15]))  # IL15Ra
                        actVec = Op(unkVec)

                        # account for pSTAT5 saturation and then normalize from 0 to 1
                        actVec = actVec / (actVec + scale)
                        actVec = actVec / T.max(actVec)

                        # concatenate the cell's IL-2 data to the IL-15 data so it matches actVec
                        data_cat = np.concatenate((self.IL2_data[(i * 4):((i + 1) * 4)], self.IL15_data[(i * 4):((i + 1) * 4)]))
                        print("data_cat[0]", data_cat[0])

                        newVec = T.reshape(actVec, data_cat.shape)  # TODO: use print statements to make sure that this reshape has been done correctly

                        temp = newVec - data_cat  # append residual to the list
                        print(temp)
                        tot_res = T.set_subtensor(tot_res[i], temp)

        return tot_res
                    

class build_model:
    """ Build the overall model handling Visterra fitting. """

    def __init__(self, traf=True):
        self.traf = traf
        self.act = pSTAT_activity()
        self.M = self.build()
        
    def build(self):
        M = pm.Model()

        with M:
            if self.traf:
                kfwd, endo, activeEndo, kRec, kDeg, sortF = commonTraf()
            else:
                kfwd = pm.Lognormal('kfwd', mu=np.log(0.001), sd=0.5, shape=1)
                # Assigning trafficking to zero to fit without trafficking
                endo = T.zeros(1, dtype=np.float64)
                activeEndo = T.zeros(1, dtype=np.float64)
                kRec = T.zeros(1, dtype=np.float64)
                kDeg = T.zeros(1, dtype=np.float64)
                sortF = T.ones(1, dtype=np.float64) * 0.5
            
            rxnrates = pm.Lognormal('rxn', sd=0.5, shape=6)  # 6 reverse rxn rates for IL2/IL15
            nullRates = T.ones(4, dtype=np.float64)  # k27rev, k31rev, k33rev, k35rev
            scale = pm.Lognormal('scales', mu=np.log(100.), sd=1, shape=1)  # create scaling constant for activity measurements
            
            # plug in 0 for all receptor rates... will override with true measurements in act.calc()
            Rexpr = T.zeros(4, dtype=np.float64)

            unkVec = T.concatenate((kfwd, rxnrates, nullRates, endo, activeEndo, sortF, kRec, kDeg, Rexpr, nullRates * 0.0))

            Y_act = self.act.calc(unkVec, scale)  # fitting the data based on dst15.calc for the given parameters
            sd_act = T.minimum(T.std(Y_act), 0.03)  # Add bounds for the stderr to help force the fitting solution
            pm.Deterministic('Y_act', T.sum(T.square(Y_act)))
            pm.Normal('fitD_act', sd=sd_act, observed=Y_act)  # experimental-derived stderr is used

            # Save likelihood
            pm.Deterministic('logp', M.logpt)

        return M
