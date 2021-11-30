"""
This file contains functions that are used in multiple figures.
"""
import sys
import logging
import time
from string import ascii_lowercase
import matplotlib
import seaborn as sns
import numpy as np
import svgutils.transform as st
from matplotlib import gridspec, pyplot as plt

matplotlib.use('AGG')
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 2
matplotlib.rcParams["ytick.major.pad"] = 2
matplotlib.rcParams["xtick.minor.pad"] = 1.9
matplotlib.rcParams["ytick.minor.pad"] = 1.9
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["legend.borderpad"] = 0.35


dosemat = np.array([84, 28, 9.333333, 3.111, 1.037037, 0.345679, 0.115226, 0.038409, 0.012803, 0.004268, 0.001423, 0.000474])


def getSetup(figsize, gridd, multz=None, empts=None):
    """ Establish figure set-up with subplots. """
    sns.set(style="whitegrid", font_scale=0.7, color_codes=True, palette="colorblind", rc={"grid.linestyle": "dotted", "axes.linewidth": 0.6})

    # create empty list if empts isn't specified
    if empts is None:
        empts = []

    if multz is None:
        multz = dict()

    # Setup plotting space and grid
    f = plt.figure(figsize=figsize, constrained_layout=True)
    gs1 = gridspec.GridSpec(*gridd, figure=f)

    # Get list of axis objects
    x = 0
    ax = list()
    while x < gridd[0] * gridd[1]:
        if x not in empts and x not in multz.keys():  # If this is just a normal subplot
            ax.append(f.add_subplot(gs1[x]))
        elif x in multz.keys():  # If this is a subplot that spans grid elements
            ax.append(f.add_subplot(gs1[x: x + multz[x] + 1]))
            x += multz[x]
        x += 1

    return (ax, f)


def subplotLabel(axs):
    """ Place subplot labels on figure. """
    for ii, ax in enumerate(axs):
        ax.text(-0.2, 1.25, ascii_lowercase[ii], transform=ax.transAxes, fontsize=16, fontweight="bold", va="top")


def overlayCartoon(figFile, cartoonFile, x, y, scalee=1, scale_x=1, scale_y=1):
    """ Add cartoon to a figure file. """

    # Overlay Figure cartoons
    template = st.fromfile(figFile)
    cartoon = st.fromfile(cartoonFile).getroot()

    cartoon.moveto(x, y)
    cartoon.scale(scalee)

    template.append(cartoon)
    template.save(figFile)


def genFigure():
    """ Main figure generation function. """
    fdir = './output/'
    start = time.time()
    nameOut = 'figure' + sys.argv[1]

    exec('from ckine.figures import ' + nameOut)
    ff = eval(nameOut + '.makeFigure()')
    ff.savefig(fdir + nameOut + '.svg', dpi=ff.dpi, bbox_inches='tight', pad_inches=0)

    if sys.argv[1] == 'C2':
        # Overlay Figure 2 cartoon
        overlayCartoon(fdir + 'figureC2.svg',
                       './ckine/graphics/tensor4D.svg', 5, 6, scalee=1.62)

    if sys.argv[1] == 'C3':
        # Overlay Figure 3 cartoon
        overlayCartoon(fdir + 'figureC3.svg',
                       './ckine/graphics/ModelCartoon.svg', 1450, 0, scalee=0.023)

    logging.info('%s is done after %s seconds.', nameOut, time.time() - start)


cellSTATlimDict = {"Treg": (0, 60000),
                   "Thelper": (0, 40000),
                   "CD8": (0, 30000),
                   "NK": (0, 8000)}


def plotDoseResponses(ax, df, mut, cellType, val=False):
    """Plots all experimental vs. Predicted Values"""
    if isinstance(cellType, str):
        if val:
            expData = df.loc[(df.Ligand == mut) & (df.Valency == val) & (df.Cell == cellType)]
        else:
            expData = df.loc[(df.Ligand == mut) & (df.Cell == cellType)]
        valList = expData.Valency.unique()

        if val:
            if val == 1:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", ax=ax, hue="Valency", palette=["darkblue"])
                sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", ax=ax, hue="Valency", palette=["darkblue"])
                ax.set(title=cellType, xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
            if val == 2:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", ax=ax, hue="Valency", palette=["springgreen"])
                sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", ax=ax, hue="Valency", palette=["springgreen"])
                ax.set(title=cellType, xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
        else:
            if len(valList) > 1:
                sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", hue="Valency", ax=ax, legend="brief", palette=["darkblue", "seagreen"])
                sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", hue="Valency", ax=ax, legend="brief", palette=["darkblue", "seagreen"])
                ax.set(title=cellType, xlabel=r"$log_{10}$ " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                handles, labels = ax.get_legend_handles_labels()
                ax.legend([handles[0]] + handles[4::], [labels[0]] + labels[4::])
            else:
                if valList[0] == 1:
                    sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", ax=ax, hue="Valency", palette=["darkblue"])
                    sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", ax=ax, color="darkblue")
                    ax.set(title=cellType, xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[0:2] + handles[4::], labels[0:2] + labels[4::])
                if valList[0] == 2:
                    sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", style="Time", ax=ax, hue="Valency", palette=["seagreen"])
                    sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", ax=ax, color="seagreen")
                    ax.set(title=cellType, xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType])
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[0:2] + handles[4::], labels[0:2] + labels[4::])
    else:
        expData = df.loc[(df.Ligand == mut) & (df.Valency == val) & (df.Cell.isin(cellType))]
        sns.scatterplot(x="Dose", y="Experimental", data=expData, label="Experimental", hue="Cell", ax=ax)
        sns.lineplot(x="Dose", y="Predicted", data=expData, label="Predicted", hue="Cell", ax=ax)
        if val == 1:
            ax.set(title=cellType[0] + "s", xlabel=r"$log_{10}$ Monomeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType[0]])
        if val == 2:
            ax.set(title=cellType[0] + "s", xlabel=r"$log_{10}$ Dimeric " + mut + " (nM)", ylabel="pSTAT", xscale="log", xlim=(1e-4, 1e2), ylim=cellSTATlimDict[cellType[0]])
