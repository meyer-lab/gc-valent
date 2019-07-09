"""
This creates Figure S2. Full panel of Geweke convergence tests.
"""
import string
import pymc3 as pm
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import seaborn as sns
from .figureCommon import subplotLabel, getSetup, traf_names
from ..imports import import_samples_2_15, import_samples_4_7


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((12, 6), (3, 4))

    # Add subplot labels
    for ii, item in enumerate(ax):
        subplotLabel(item, string.ascii_uppercase[ii])

    plot_geweke_2_15(ax[0], True)
    plot_geweke_2_15(ax[4], False)
    plot_geweke_4_7(ax[8:11])

    return f


def plot_geweke_2_15(ax, traf):
    """ Uses geweke criterion to evaluate model convergence during fitting. """
    trace = import_samples_2_15(Traf=traf, ret_trace=True)  # return the trace

    # use use trace to calculate geweke z-scores
    score = pm.diagnostics.geweke(trace, first=0.1, last=0.5, intervals=20)

    # take score from the 10th interval in the chain... can change this to a random int later
    dictt = {r'$k_{4}$': score[0]['rxn'][0][10, 1], 
             r'$k_{5}$': score[0]['rxn'][1][10, 1],
             r'$k_{16}$': score[0]['rxn'][2][10, 1],
             r'$k_{17}$': score[0]['rxn'][3][10, 1],
             r'$k_{22}$': score[0]['rxn'][4][10, 1],
             r'$k_{23}$': score[0]['rxn'][5][10, 1],
             'IL-2Rα': score[0]['Rexpr_2Ra'][10, 1],
             'IL-2Rβ': score[0]['Rexpr_2Rb'][10, 1],
             'IL-15Rα': score[0]['Rexpr_15Ra'][10, 1],
             r'$C_{5}$': score[0]['scales'][10, 1],
             r'$k_{fwd}$': score[0]['kfwd'][10, 1]}

    if traf:  # add the trafficking parameters if necessary & set proper title
        dictt.update({r'$k_{endo}$': score[0]['endo'][10, 1],
                     r'$k_{endo,a}$': score[0]['activeEndo'][10, 1],
                     r'$f_{sort}$': score[0]['sortF'][10, 1],
                     r'$k_{rec}$': score[0]['kRec'][10, 1],
                     r'$k_{deg}$': score[0]['kDeg'][10, 1]})
        ax.set_title(r'IL-2/-15 trafficking model')
    else:
        ax.set_title(r'IL-2/-15 no trafficking model')

    df = pd.DataFrame.from_dict(dictt, orient='index')
    sns.scatterplot(data=np.abs(df), ax=ax)

    ax.set_xticklabels(list(dictt.keys()),  # use keys from dict as x-axis labels
                       rotation=25,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    ax.get_legend().set_visible(False)  # remove legend created by sns
    ax.axhline(1., c='r')  # line to denote acceptable threshold of standard deviations
    ax.set(ylim=(-0.1, 1.25))



def plot_geweke_4_7(ax, traf=True):
    """ Generating Geweke plots using the traces from IL-4 and IL-7 fitting to Gonnord data. """
    trace = import_samples_4_7(ret_trace=True)  # return the trace

    # use use trace to calculate geweke z-scores
    score = pm.diagnostics.geweke(trace, first=0.1, last=0.5, intervals=20)

    # take score from the 10th interval in the chain... can change this to a random int later
    dictt = {r'$k_{27}$': score[0]['k27rev'][10, 1], 
             r'$k_{33}$': score[0]['k33rev'][10, 1],
             r'$C_{5}$': score[0]['scales'][0][10, 1],
             r'$C_{6}$': score[0]['scales'][1][10, 1],
             r'$k_{fwd}$': score[0]['kfwd'][10, 1]}

    if traf:  # add the trafficking parameters if necessary & set proper title
        dictt.update({r'$k_{endo}$': score[0]['endo'][10, 1],
                     r'$k_{endo,a}$': score[0]['activeEndo'][10, 1],
                     r'$f_{sort}$': score[0]['sortF'][10, 1],
                     r'$k_{rec}$': score[0]['kRec'][10, 1],
                     r'$k_{deg}$': score[0]['kDeg'][10, 1]})
        ax.set_title(r'IL-4/-7 trafficking model')
    else:
        ax.set_title(r'IL-4/-7 no trafficking model')

    df = pd.DataFrame.from_dict(dictt, orient='index')
    sns.scatterplot(data=np.abs(df), ax=ax)

    ax.set_xticklabels(list(dictt.keys()),  # use keys from dict as x-axis labels
                       rotation=25,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    ax.get_legend().set_visible(False)  # remove legend created by sns
    ax.axhline(1., c='r')  # line to denote acceptable threshold of standard deviations
    ax.set(ylim=(-0.1, 1.25))

