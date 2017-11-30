from model import dy_dt_IL15_wrapper
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pds
from theano.compile.ops import as_op
import theano.tensor as T
import pymc3 as pm
import pickle as pk
import bz2
import concurrent.futures

pool = concurrent.futures.ProcessPoolExecutor()

# this just takes the output of odeint (y values) and determines pSTAT activity
def IL15_pSTAT_activity(ys):
    # pSTAT activation is based on sum of IL15_IL2Rb_gc and IL15_IL15Ra_IL2Rb_gc  and IL15-IL15Ra-IL2Rb and IL15-IL15Ra-gc at long time points
    activity = ys[1,14] + ys[1,15] + ys[1,16] + ys[1,17]
    return activity

# this takes the values of input parameters and calls odeint, then puts the odeint output into IL15_pSTAT_activity
def IL15_activity_input(y0, t, IL15, k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev):
    args = (IL15, k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev)
    ts = np.linspace(0., t, 2)
    ys = odeint(dy_dt_IL15_wrapper, y0, ts, args, mxstep = 6000)
    act = IL15_pSTAT_activity(ys)
    return act


# this takes all the desired IL15 values we want to test and gives us the maximum activity value
# IL15 values pretty much ranged from 5 x 10**-4 to 500 nm with 8 points in between
@as_op(itypes=[T.dscalar, T.dscalar, T.dscalar, T.dscalar, T.dscalar, T.scalar], otypes=[T.dmatrix])

def IL15_activity_values(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev):    
    y0 = np.array([1000.,1000.,1000., 0., 0., 0., 0., 0., 0., 0.])

    t = 50.
    IL15s = np.logspace(-3.3, 2.7, 8) # 8 log-spaced values between our two endpoints
    table = np.zeros((8, 2))
    output = list()

    for ii in range(len(IL15s)):
        output.append(pool.submit(IL15_activity_input, y0, t, IL15s[ii], k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev))

    for ii in range(len(IL15s)):
        table[ii, 1] = output[ii].result()
    
    table[:, 0] = IL15s

    return table


def IL15_percent_activity(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev):
    values = IL15_activity_values(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev)
    maximum = T.max(values[:,1], 0) # find the max value in the second column for all rows

    new_table = T.stack((values[:, 0], 100. * values[:, 1] / maximum), axis=1) # IL15 values in first column are the same
    # activity values in second column are converted to percents relative to maximum
    
    return new_table


def plot_IL15_percent_activity(y0, t, k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev):
    new_table = IL15_percent_activity(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev)

    x = math.log10(new_table[:, 0]) # changing the x values to the log10(nM) values that were in the published graph

    plt.rcParams.update({'font.size': 8})
    plt.xlabel("IL15 concentration (log(nm))")
    plt.ylabel("percent activation of pSTAT")
    plt.scatter(x[:], new_table[:,1])
    plt.show()

class IL15_sum_squared_dist:
    
    def load(self):
        data = pds.read_csv("./data/IL2_IL15_extracted_data.csv") # imports csv file into pandas array
        self.numpy_data = data.as_matrix() #the IL2_IL2Ra- data is within the 3rd column (index 2)
        
    def calc(self, k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev):
        activity_table = IL15_percent_activity(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev)
        diff_data = self.numpy_data[:,7] - activity_table[:,1]
        return np.squeeze(diff_data)
    
def store_data(class_name, fit_results):
    x = pk.dump(class_name, bz2.BZ2File(fit_results + '.pkl', 'wb'))
    return x

class build_model:
    
    def __init__(self):
        self.dst = IL15_sum_squared_dist()
        self.dst.load()
        
    def build(self):   
        self.M = pm.Model()
        
        with self.M:
            k13fwd = pm.Lognormal('k13fwd', mu=0, sd=3) # we need to add a standard deviation and they're all based on a lognormal scale
            k15rev = pm.Lognormal('k15rev', mu=0, sd=3)
            k17rev = pm.Lognormal('k17rev', mu=0, sd=3)
            k18rev = pm.Lognormal('k18rev', mu=0, sd=3)
            k22rev = pm.Lognormal('k22rev', mu=0, sd=3)
            k23rev = pm.Lognormal('k23rev', mu=0, sd=3)
            
            Y = self.dst.calc(k13fwd, k15rev, k17rev, k18rev, k22rev, k23rev)
            
            Y_obs = pm.Normal('fitD', mu=0, sd=T.std(Y), observed=Y)
        
        return Y_obs
    
    def sampling(self):
        with self.M:
            start = pm.find_MAP()
            step = pm.Metropolis()
            self.trace = pm.sample(5000, step, start=start) # original value should be 5 to shorten time
            
            
M = build_model()
M.build()
M.sampling()

# _ = plt.hist(build_model.M.trace['k4fwd'],100)

store_data(M, "IL15_model_results")