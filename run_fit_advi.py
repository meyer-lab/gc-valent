"""
This is a pseudo-replica of run_fit.py. This calls fit_advi.py instead of fit.py.
"""

# Set matplotlib backend so python remains in the background
import matplotlib
matplotlib.use("Agg")
from ckine.fit_advi import build_model
import pymc3 as pm

if __name__ == "__main__": #only go into this loop if you're running fit_advi.py directly instead of running a file that calls fit.py
    M = build_model()
    M.build()
    M.sampling()
    pm.backends.text.dump("IL2_IL15_advi_results", M.trace) #instead of pickling data we dump it into file that can be accessed by read_fit_data.py
    