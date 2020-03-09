"""
This file contains functions that are used in multiple figures.
"""
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.cm as cm
from matplotlib import gridspec, pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from ..imports import import_pstat
from ..flow import nllsq, exp_dec


matplotlib.rcParams['legend.labelspacing'] = 0.2
matplotlib.rcParams['legend.fontsize'] = 8
matplotlib.rcParams['xtick.major.pad'] = 2
matplotlib.rcParams['ytick.major.pad'] = 2
matplotlib.rcParams['xtick.minor.pad'] = 1.9
matplotlib.rcParams['ytick.minor.pad'] = 1.9
matplotlib.rcParams['legend.handletextpad'] = 0.5
matplotlib.rcParams['legend.handlelength'] = 0.5
matplotlib.rcParams['legend.framealpha'] = 0.5
matplotlib.rcParams['legend.markerscale'] = 0.7
matplotlib.rcParams['legend.borderpad'] = 0.35


def getSetup(figsize, gridd, multz=None, empts=None):
    """ Establish figure set-up with subplots. """
    sns.set(style="whitegrid",
            font_scale=0.7,
            color_codes=True,
            palette="colorblind",
            rc={'grid.linestyle': 'dotted',
                'axes.linewidth': 0.6})

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
            ax.append(f.add_subplot(gs1[x:x + multz[x] + 1]))
            x += multz[x]
        x += 1

    return (ax, f)


def set_bounds(ax):
    """ Set bounds of component plots. """
    x_max = np.max(np.absolute(np.asarray(ax.get_xlim()))) * 1.1
    y_max = np.max(np.absolute(np.asarray(ax.get_ylim()))) * 1.1

    ax.set_xlim(-x_max, x_max)
    ax.set_ylim(-y_max, y_max)


def subplotLabel(ax, letter, hstretch=1, ystretch=1):
    """ Label each subplot """
    ax.text(-0.2 / hstretch, 1.2 / ystretch, letter, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top')


def traf_names():
    """ Returns a list of the trafficking parameters in order they appear within unkVec. """
    return [r'$\mathrm{k_{endo}}$', r'$\mathrm{k_{endo,a}}$', r'$\mathrm{k_{rec}}$', r'$\mathrm{k_{deg}}$']


def plot_conf_int(ax, x_axis, y_axis, color, label=None):
    """ Shades the 25-75 percentiles dark and the 10-90 percentiles light. The percentiles are found along axis=1. """
    y_axis_top = np.percentile(y_axis, 90., axis=1)
    y_axis_bot = np.percentile(y_axis, 10., axis=1)
    ax.fill_between(x_axis, y_axis_top, y_axis_bot, color=color, alpha=0.4)

    y_axis_top = np.percentile(y_axis, 75., axis=1)
    y_axis_bot = np.percentile(y_axis, 25., axis=1)
    ax.fill_between(x_axis, y_axis_top, y_axis_bot, color=color, alpha=0.7, label=label)
    if label is not None:
        ax.legend()


def plot_cells(ax, factors, component_x, component_y, cell_names):
    """This function plots the combination decomposition based on cell type."""
    colors = cm.rainbow(np.linspace(0, 1, len(cell_names)))
    markersCells = ['^', '*', 'D', 's', 'X', 'o', '4', 'H', 'P', '*', 'D', 's', 'X']

    for ii, _ in enumerate(factors[:, component_x - 1]):
        ax.scatter(factors[ii, component_x - 1], factors[ii, component_y - 1], c=[colors[ii]], marker=markersCells[ii], label=cell_names[ii])

    ax.set_title('Cells')
    ax.set_xlabel('Component ' + str(component_x))
    ax.set_ylabel('Component ' + str(component_y))
    ax.set_xlim(left=-0.03)
    ax.legend()


def overlayCartoon(figFile, cartoonFile, x, y, scalee=1, scale_x=1, scale_y=1):
    """ Add cartoon to a figure file. """
    import svgutils.transform as st

    # Overlay Figure cartoons
    template = st.fromfile(figFile)
    cartoon = st.fromfile(cartoonFile).getroot()

    cartoon.moveto(x, y, scale=scalee)
    cartoon.scale_xy(scale_x, scale_y)

    template.append(cartoon)
    template.save(figFile)


def plot_ligands(ax, factors, ligand_names, cutoff=0.0):
    """Function to put all ligand decomposition plots in one figure."""
    ILs, _, _, _, _ = import_pstat()  # Cytokine stimulation concentrations in nM
    n_ligands = len(ligand_names)
    ILs = np.flip(ILs)
    colors = sns.color_palette()
    legend_shape = []
    markers = ['.', '^', 'd', '*']

    for ii, name in enumerate(ligand_names):
        legend_shape.append(Line2D([0], [0], color='k', marker=markers[ii], label=name, linestyle=''))  # Make ligand legend elements

    for ii in range(factors.shape[1]):
        componentLabel = True
        for jj in range(n_ligands):
            idx = range(jj * len(ILs), (jj + 1) * len(ILs))

            # If the component value never gets over cutoff, then don't plot the line
            if np.max(factors[idx, ii]) > cutoff:
                if componentLabel:
                    ax.plot(ILs, factors[idx, ii], color=colors[ii], label='Cmp. ' + str(ii + 1))
                    componentLabel = False
                else:
                    ax.plot(ILs, factors[idx, ii], color=colors[ii])
                ax.scatter(ILs, factors[idx, ii], color=colors[ii], marker=markers[jj])

    ax.add_artist(ax.legend(handles=legend_shape, loc=2))

    ax.set_xlabel('Ligand Concentration (nM)')
    ax.set_ylabel('Component')
    ax.set_xscale('log')
    ax.set_title('Ligands')

    # Place legend
    ax.legend(loc=6)


def plot_timepoints(ax, ts, factors):
    """Function to put all timepoint plots in one figure."""
    colors = sns.color_palette()
    for ii in range(factors.shape[1]):
        ax.plot(ts, factors[:, ii], c=colors[ii], label='Component ' + str(ii + 1))

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Component')
    ax.set_title('Time')
    ax.legend()


def legend_2_15(ax, location="center right"):
    """ Plots a legend for all the IL-2 and IL-15 related plots in its own subpanel. """
    legend_elements = [Patch(facecolor='darkorchid', label='IL-2'),
                       Patch(facecolor='goldenrod', label='IL-15'),
                       Line2D([0], [0], marker='o', color='w', label='IL-2Rα+',
                              markerfacecolor='k', markersize=16),
                       Line2D([0], [0], marker='^', color='w', label='IL-2Rα-',
                              markerfacecolor='k', markersize=16)]
    ax.legend(handles=legend_elements, loc=location, prop={'size': 16})
    ax.axis('off')  # remove the grid


def plot_hist(axes, sample, channels):
    """ Plots histogram of signal for each well/channel in a sample. """
    for i, s in enumerate(sample):
        tform = s.transform('hlog', channels=channels[i])
        data = tform.data[[channels[i]]][0:]
        axes[i].hist(data[channels[i]], bins=100)


def plot_fsc_ssc(axes, sample):
    """ Plots forward and side scatter for a given sample (all wells). """
    for i, s in enumerate(sample):
        s.plot(['SSC-H', 'FSC-H'], ax=axes[i])
        x0, x1 = axes[i].get_xlim()
        y0, y1 = axes[i].get_ylim()
        axes[i].set_aspect((x1 - x0) / (y1 - y0))


def plot_regression(ax, sample, channels, receptors, recQuant, first=0, skip=False):
    """ Plots regression of signal to bead capacity. """
    means = np.zeros(len(recQuant))
    for i, s in enumerate(sample):
        if skip:
            if i < first:
                continue
        tform = s.transform('hlog', channels=channels[i - first])
        data = tform.data[[channels[i - first]]][0:]
        avg_signal = np.mean(data[str(channels[i - first])])
        means[i - first] = avg_signal
    ax.scatter(recQuant, means)
    lsq = nllsq([1, 1000], recQuant, means)
    xs = np.linspace(np.amin(recQuant), np.amax(recQuant), num=1000)
    ax.plot(xs, exp_dec(lsq, xs))
    ax.set_xlabel('Bead Capacity')
    ax.set_ylabel('Average Signal (' + str(receptors[4 + first]) + ')')
