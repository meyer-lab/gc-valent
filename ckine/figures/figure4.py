"""
This creates Figure 4.
"""
import string
import numpy as np
import math
from scipy.optimize import least_squares
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 4))

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
    
    ckineConc, cell_names_pstat, IL2_data, _ = import_pstat()
    
    x0 = [10**-10, 2.]
    
    for i, _ in enumerate(cell_names_pstat):
        test = IL2_data[(i * 4):((i + 1) * 4)]
        nllsq(ax[i], x0, ckineConc.astype(np.float), test[3])

    return f


def nllsq(ax, x0, xdata, ydata):
    lsq_res = least_squares(residuals, x0, args = (xdata, ydata), bounds=([10.**-13., 0.],[10.**-7, 5.]), verbose=2)
    print(lsq_res.x)
    y = hill_equation(xdata, lsq_res.x)
    ax.scatter(xdata, ydata)
    ax.plot(xdata, y)

def hill_equation(x, x0):
    k = x0[0]
    n = x0[1]
    y = (x**n)/((k**n) + (x**n))
    return y

def residuals(x0, x, y):
    res = y - hill_equation(x, x0)
    return res

