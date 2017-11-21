from model import dy_dt_IL2_wrapper
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pds
from theano.compile.ops import as_op
import theano.tensor as T
from pymc3 import Model, Lognormal, Normal

# this just takes the output of odeint (y values) and determines pSTAT activity
def IL2_pSTAT_activity(ys):
    # pSTAT activation is based on sum of IL2_IL2Rb_gc and IL2_IL2Ra_IL2Rb_gc at long time points
    activity = ys[1, 8] + ys[1, 9]
    return activity

# this takes the values of input parameters and calls odeint, then puts the odeint output into IL2_pSTAT_activity
def IL2_activity_input(y0, t, IL2, k4fwd, k5rev, k6rev):
    args = (IL2, k4fwd, k5rev, k6rev)
    ts = np.linspace(0., t, 2)
    ys = odeint(dy_dt_IL2_wrapper, y0, ts, args, mxstep = 6000)
    act = IL2_pSTAT_activity(ys)
    return act

# this takes all the desired IL2 values we want to test and gives us the maximum activity value
# IL2 values pretty much ranged from 5 x 10**-4 to 500 nm with 8 points in between
def IL2_activity_values(y0, t, k4fwd, k5rev, k6rev):
    IL2s = np.logspace(-3.3, 2.7, 8) # 8 log-spaced values between our two endpoints
    activity = np.zeros(8)
    table = np.zeros((8,2))
    for ii in range (IL2s.shape[0]):
        activity[ii] = IL2_activity_input(y0, t, IL2s[ii], k4fwd, k5rev, k6rev)
    
    table[:, 0] = IL2s
    table[:, 1] = activity

    return table


def IL2_percent_activity(y0, t, k4fwd, k5rev, k6rev):
    values = IL2_activity_values(y0, t, k4fwd, k5rev, k6rev)
    maximum = np.amax(values[:,1], 0) # find the max value in the second column for all rows
    new_table = np.zeros((8,2))

    new_table[:, 0] = values[:, 0] # IL2 values in first column are the same
    new_table[:, 1] = 100. * values[:, 1] / maximum # activity values in second column are converted to percents relative to maximum
    
    return new_table


def plot_IL2_percent_activity(y0, t, k4fwd, k5rev, k6rev):
    new_table = IL2_percent_activity(y0, t, k4fwd, k5rev, k6rev)

    x = math.log10(new_table[:, 0]) # changing the x values to the log10(nM) values that were in the published graph

    plt.rcParams.update({'font.size': 8})
    plt.xlabel("IL2 concentration (log(nm))")
    plt.ylabel("percent activation of pSTAT")
    plt.scatter(x[:], new_table[:,1])
    plt.show()


@as_op(itypes=[T.dscalar, T.dscalar, T.dscalar], otypes=[T.dvector])
def IL2_sum_squared_distance(k4fwd, k5rev, k6rev):
    y0 = np.array([1000.,1000.,1000.,0.,0.,0.,0.,0.,0.,0.])
    t = 50.

    activity_table = IL2_percent_activity(y0, t, k4fwd, k5rev, k6rev) # generates output from percent activity function
    data = pds.read_csv("./data/IL2_IL15_extracted_data.csv") # imports csv file into pandas array
    
    # trying to perform calculation using numpy arrays
    numpy_data = data.as_matrix() #the IL2_IL2Ra- data is within the 3rd column (index 2)

    diff_data = numpy_data[:,6] - activity_table[:,1] # second column represents IL2_IL2Ra+ data
    # deciding to return diff_data for now

    return np.squeeze(diff_data)

#print (IL2_sum_squared_distance([1000.,1000.,1000.,0.,0.,0.,0.,0.,0.,0.], 50., 1., 1., 1.))


with Model():
    k4fwd = Lognormal('k4fwd', mu=0, sd=3) # do we need to add a standard deviation? Yes, and they're all based on a lognormal scale
    k5rev = Lognormal('k5rev', mu=0, sd=3)
    k6rev = Lognormal('k6rev', mu=0, sd=3)
    
    Y = IL2_sum_squared_distance(k4fwd, k5rev, k6rev)
    
    Y_obs = Normal('fitD', mu=0, sd=T.std(Y), observed=Y)
    