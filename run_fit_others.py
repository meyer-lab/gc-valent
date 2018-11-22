#!/usr/bin/env python3

# Set matplotlib backend so python remains in the background
import matplotlib
import numpy as np
import random
matplotlib.use("Agg")
from ckine.fit_others import build_model
import pymc3 as pm

if __name__ == "__main__": #only go into this loop if you're running fit.py directly instead of running a file that calls fit.py
    M = build_model(pretreat=True)
    M.build()
    rates = np.logspace(-3, 3)
    scales = np.logspace(-2, 5)
    # randomly select arguments for sampling from ex
    for ii in range(25):
        k27 = random.choice(rates)
        k33 = random.choice(rates)
        scale1 = random.choice(scales)
        scale2 = random.choice(scales)
        print("start - k27rev: " + str(k27) + ", k33rev: " + str(k33) + ", scale1: " + str(scale1) + ", scale2: " + str(scale2))
        M.sampling(k27, k33, scale1, scale2)
        print(M.trace)
    # pm.backends.text.dump("IL4-7_model_results", M.trace) #instead of pickling data we dump it into file that can be accessed by read_fit_data.py
    # pm.backends.text.Text("IL4-7_find_MAP_results", model=M)
