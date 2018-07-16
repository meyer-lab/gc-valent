"""

"""
import seaborn as sns
from matplotlib import gridspec, pyplot as plt

def getSetup(figsize, gridd):
    sns.set(style="whitegrid",
            font_scale=0.7,
            color_codes=True,
            palette="colorblind",
            rc={'grid.linestyle':'dotted',
                'axes.linewidth':0.6})

    # Setup plotting space
    f = plt.figure(figsize=figsize)

    # Make grid
    gs1 = gridspec.GridSpec(*gridd)

    # Get list of axis objects
    ax = [f.add_subplot(gs1[x]) for x in range(gridd[0] * gridd[1])]

    return (ax, f)


def subplotLabel(ax, letter, hstretch=1):
    ax.text(-0.2 / hstretch, 1.2, letter, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top')

def rev_rxn_names():
    return ['k4rev', 'k5rev', 'k16rev', 'k17rev', 'k22rev', 'k23rev', 'k27rev', 'k31rev', 'k33rev', 'k35rev']

def traf_names():
    return ['endo', 'activeEndo', 'sortF', 'kRec', 'kDeg']

def Rexpr_names():
    return ['IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra', 'IL21Ra']
