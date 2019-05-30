"""
This creates Figure 4.
"""
import string
import numpy as np
import math
from scipy.optimize import curve_fit
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
    
    for i, _ in enumerate(cell_names_pstat):
        test = IL2_data[(i * 4):((i + 1) * 4)]
        nllsq(ax[i], ckineConc.astype(np.float), test[3])

    return f


def nllsq(ax, xdata, ydata):
    popt, pcov = curve_fit(hill_equation, xdata, ydata, bounds=([10.**-13., 0.],[10.**-7, 1.]))
    y = hill_equation(xdata, *popt)
    ax.scatter(xdata, ydata)
    ax.plot(xdata, y)

def hill_equation(x, k, n):
    y = (x**n)/((k**n) + (x**n))
    return y

"""
def residuals(x, y, k, n):
    res = y - hill_equation(x, k, n)
    return res
"""
