"""
This creates Figure 4. CP decomposition of measured pSTAT data.
"""
import string
import numpy as np
import math
from scipy.optimize import least_squares, minimize
from .figureCommon import subplotLabel, getSetup
from ..imports import import_pstat


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (4, 4))

    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])
    
    ckineConc, cell_names_pstat, IL2_data, _ = import_pstat()
    
    x0 = [1, 2., 1000.]
    
    print(np.log10(ckineConc.astype(np.float)*10**4))
    
    
    for i, _ in enumerate(cell_names_pstat):
        test = IL2_data[(i * 4):((i + 1) * 4)]
        test2 = test[3] # add log?
        nllsq(ax[i], x0, np.log10(ckineConc.astype(np.float)*10**4), test2)

    return f


def nllsq(ax, x0, xdata, ydata):
    lsq_res = least_squares(residuals, x0, args=(xdata, ydata), bounds=([0., 0., 0.], [10., 10., 10**5.]), jac='3-point')
    y = hill_equation(xdata, lsq_res.x)
    ax.scatter(xdata, ydata)
    ax.plot(xdata, y)

def hill_equation(x, x0):
    k = x0[0]
    n = x0[1]
    A = x0[2]
    xk = np.power(x / k, n)
    return A * xk / (1.0 + xk)

def residuals(x0, x, y):
    return hill_equation(x, x0) - y
