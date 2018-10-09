"""
This creates Figure 4.
"""
import string
from .figureCommon import subplotLabel, getSetup
from ..fit_others import build_model as build_model_others
from ..fit import build_model


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    f.tight_layout()

    return f


def relativeGC(ax, unkVec):
    """ This function compares the relative complex affinities for GC. The rates included in this violing plot will be k4rev, k10rev, k17rev, k22rev, k27rev, and k33rev. We're currently ignoring k31rev (IL9) and k35rev (IL21) since we don't fit to any of its data. """

    # assign values from unkVec
    k4rev, k5rev, k16rev, k17rev, k22rev, k23rev, k27rev, k33rev = unkVec[7, :], unkVec[8, :], unkVec[9, :], unkVec[10, :], unkVec[11, :], unkVec[12, :], unkVec[13, :], unkVec[15, :]
    
    # back-out k10 with ratio
    k10rev = 12.0 * k5rev / 1.5 # doi:10.1016/j.jmb.2004.04.038
    
    # add each rate duo as separate column in dataframe
    df = pd.DataFrame({'k4': k4rev, 'k10': k10rev, 'k17': k17rev, 'k22': k22rev, 'k27': k27rev, 'k33': k33rev})
    
    # df.columns = ['IL2Ra', 'IL2Rb', 'gc', 'IL15Ra']
    a = sns.violinplot(data=np.log10(df), ax=ax, linewidth=0, bw=10)
    a.set_xticklabels(a.get_xticklabels(), rotation=40, rotation_mode="anchor", ha="right", fontsize=8, position=(0, 0.075))
    a.set(title="Relative gc affinity", ylabel="log10 of value")
    