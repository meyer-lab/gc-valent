#!/usr/bin/env python3

# Set matplotlib backend so python remains in the background
import matplotlib
matplotlib.use("Agg")
from ckine.fit_others import build_model
import pymc3 as pm

if __name__ == "__main__": #only go into this loop if you're running fit.py directly instead of running a file that calls fit.py
    M = build_model(pretreat=True)
    M.build()
    ex = np.linspace(10**-8, 10**4)
    # randomly select arguments for sampling from ex
    M.sampling()
    print(M.trace)
    # pm.backends.text.dump("IL4-7_model_results", M.trace) #instead of pickling data we dump it into file that can be accessed by read_fit_data.py
    # pm.backends.text.Text("IL4-7_find_MAP_results", model=M)
