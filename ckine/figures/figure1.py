"""
This creates Figure 1.
"""
from .figureCommon import subplotLabel, getSetup, traf_names
from ..plot_model_prediction import surf_IL2Rb, pstat, surf_gc
from ..model import nParams
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import pymc3 as pm, os
from os.path import join
from ..fit import build_model



def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((10, 9), (4, 4), mults=[2], multz={2: 2})

    subplotLabel(ax[0], 'A')
    unkVec = import_samples()
    #pstat_act(ax[0:2], unkVec)
    #rateComp(ax[2], unkVec)
    surf_perc(ax[3:7], 'IL2Rb', unkVec)
    surf_perc(ax[7:11], 'gc', unkVec)
    #violinPlots(ax[11:13], unkVec)

    f.tight_layout()

    return f

def plot_structure(IL2vec, IL15vec, title, ax, x_axis, data_type, species='NONE'):
    """ This function performs the plotting of pSTAT activity and surface receptor calculations given data for both IL2 and IL15 stimulation. """
    ax.set_title(title)
    ax.set_ylim(0,120)
    if (data_type=='surf'):
        ax.plot(x_axis, IL2vec, color='darkorchid', alpha=0.7)
        ax.plot(x_axis, IL15vec, color='goldenrod', alpha=0.7)
        ax.set_ylabel("Surface " + str(species) + " (% x 100)")
        ax.set_xlabel("Time (min)")
    elif (data_type=='act'):
        ax.plot(np.log10(x_axis), IL2vec, color='darkorchid', alpha=0.5)
        ax.plot(np.log10(x_axis), IL15vec, color='goldenrod', alpha=0.5)
        ax.set_ylabel('Maximal p-STAT5 (% x 100)')
        ax.set_xlabel('log10 of cytokine concentration (nM)')
    else:
        print('invalid data_type')   

def surf_perc(ax, species, unkVec):
    """ Calculates the percent of IL2Rb or gc on the cell surface over the course of 90 mins. Cell environments match those of surface IL2Rb data collected by Ring et al. """
    if (species == 'IL2Rb'):
        surf = surf_IL2Rb()
    elif (species == 'gc'):
        surf = surf_gc()
    else:
        print('not a valid species')
        return -1

    y_max = 100.
    ts = np.array([0., 2., 5., 15., 30., 60., 90.])
    size = len(ts)
    IL2_1_plus = np.zeros((size, 500))
    IL2_500_plus = IL2_1_plus.copy()
    IL2_1_minus = IL2_1_plus.copy()
    IL2_500_minus = IL2_1_plus.copy()
    IL15_1_plus = IL2_1_plus.copy()
    IL15_500_plus = IL2_1_plus.copy()
    IL15_1_minus = IL2_1_plus.copy()
    IL15_500_minus = IL2_1_plus.copy()

    for ii in range(0,500):
        output = surf.calc(unkVec[:, ii], ts) * y_max
        IL2_1_plus[:, ii] = output[0:(size)]
        IL2_500_plus[:, ii] = output[(size):(size*2)]
        IL2_1_minus[:, ii] = output[(size*2):(size*3)]
        IL2_500_minus[:, ii] = output[(size*3):(size*4)]
        IL15_1_plus[:, ii] = output[(size*4):(size*5)]
        IL15_500_plus[:, ii] = output[(size*5):(size*6)]
        IL15_1_minus[:, ii] = output[(size*6):(size*7)]
        IL15_500_minus[:, ii] = output[(size*7):(size*8)]

        #plot_structure(IL2_1_minus, IL15_1_minus, '1 nM and IL2Ra-', ax[0], ts, 'surf', species)
        #plot_structure(IL2_500_minus, IL15_500_minus, "500 nM and IL2Ra-", ax[1], ts, 'surf', species)
        #plot_structure(IL2_1_plus, IL15_1_plus, "1 nM and IL2Ra+", ax[2], ts, 'surf', species)
        #plot_structure(IL2_500_plus, IL15_500_plus, "500 nM and IL2Ra+", ax[3], ts, 'surf', species)
        
    # calculate 97.5th and 2.5th percentiles for each time pt. 
    IL2_1_plus_top = np.percentile(IL2_1_plus, 97.5, axis=1)
    IL2_1_plus_bot = np.percentile(IL2_1_plus, 2.5, axis=1)
    IL2_500_plus_top = np.percentile(IL2_500_plus, 97.5, axis=1)
    IL2_500_plus_bot = np.percentile(IL2_500_plus, 2.5, axis=1)
    IL2_1_minus_top = np.percentile(IL2_1_minus, 97.5, axis=1)
    IL2_1_minus_bot = np.percentile(IL2_1_minus, 2.5, axis=1)
    IL2_500_minus_top = np.percentile(IL2_500_minus, 97.5, axis=1)
    IL2_500_minus_bot = np.percentile(IL2_500_minus, 2.5, axis=1)
    IL15_1_plus_top = np.percentile(IL15_1_plus, 97.5, axis=1)
    IL15_1_plus_bot = np.percentile(IL15_1_plus, 2.5, axis=1)
    IL15_500_plus_top = np.percentile(IL15_500_plus, 97.5, axis=1)
    IL15_500_plus_bot = np.percentile(IL15_500_plus, 2.5, axis=1)
    IL15_1_minus_top = np.percentile(IL15_1_minus, 97.5, axis=1)
    IL15_1_minus_bot = np.percentile(IL15_1_minus, 2.5, axis=1)
    IL15_500_minus_top = np.percentile(IL15_500_minus, 97.5, axis=1)
    IL15_500_minus_bot = np.percentile(IL15_500_minus, 2.5, axis=1)
    
    # plot shaded regions representing 95% confidence interval
    ax[0].fill_between(ts, IL2_1_minus_top, IL2_1_minus_bot, color='darkorchid', alpha=0.5)
    ax[0].fill_between(ts, IL15_1_minus_top, IL15_1_minus_bot, color='goldenrod', alpha=0.5)
    ax[1].fill_between(ts, IL2_500_minus_top, IL2_500_minus_bot, color='darkorchid', alpha=0.5)
    ax[1].fill_between(ts, IL15_500_minus_top, IL15_500_minus_bot, color='goldenrod', alpha=0.5)
    ax[2].fill_between(ts, IL2_1_plus_top, IL2_1_plus_bot, color='darkorchid', alpha=0.5)
    ax[2].fill_between(ts, IL15_1_plus_top, IL15_1_plus_top, color='goldenrod', alpha=0.5)
    ax[3].fill_between(ts, IL2_500_plus_top, IL2_500_plus_bot, color='darkorchid', alpha=0.5)
    ax[3].fill_between(ts, IL15_500_plus_top, IL2_500_minus_bot, color='goldenrod', alpha=0.5)
    
    # set labels and titles
    ax[0].set()
    
        
    if (species == 'IL2Rb'):
        path = os.path.dirname(os.path.abspath(__file__))
        data_minus = pd.read_csv(join(path, "../data/IL2Ra-_surface_IL2RB_datasets.csv")).values # imports file into pandas array
        data_plus = pd.read_csv(join(path, "../data/IL2Ra+_surface_IL2RB_datasets.csv")).values # imports file into pandas array
        ax[0].scatter(data_minus[:,0], data_minus[:,1] * 10., color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # 1nM of IL2 in 2Ra-
        ax[0].scatter(data_minus[:,0], data_minus[:,2] * 10., color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # 1nM of IL15 in 2Ra-
        ax[1].scatter(data_minus[:,0], data_minus[:,5] * 10., color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # 500nM of IL2 in 2Ra-
        ax[1].scatter(data_minus[:,0], data_minus[:,6] * 10., color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # 500nM of IL15 in 2Ra-
        ax[2].scatter(data_plus[:,0], data_plus[:,1] * 10., color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # 1nM of IL2 in 2Ra+
        ax[2].scatter(data_plus[:,0], data_plus[:,2] * 10., color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # 1nM of IL15 in 2Ra+
        ax[3].scatter(data_plus[:,0], data_plus[:,5] * 10., color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # 500nM of IL2 in 2Ra+
        ax[3].scatter(data_plus[:,0], data_plus[:,6] * 10., color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # 500nM of IL15 in 2Ra+
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        ax[3].legend()

    
def pstat_act(ax, unkVec):
    """ This function generates the pSTAT activation levels for each combination of parameters in unkVec. The results are plotted and then overlayed with the values measured by Ring et al. """
    pstat5 = pstat()
    PTS = 30
    cytokC = np.logspace(-3.3, 2.7, PTS)
    y_max = 100.
    IL2_plus = np.zeros((PTS, 500))
    IL15_minus = IL2_plus.copy()
    IL15_plus = IL2_plus.copy()
    IL2_minus = IL2_plus.copy()
    
    # calculate activity for each unkVec for all conc.
    for ii in range(0,500):
        output = pstat5.calc(unkVec[:, ii], cytokC) * y_max
        IL2_plus[:, ii] = output[0:PTS]
        IL2_minus[:, ii] = output[PTS:(PTS*2)]
        IL15_plus[:, ii] = output[(PTS*2):(PTS*3)]
        IL15_minus[:, ii] = output[(PTS*3):(PTS*4)]

    # calculate 97.5th and 2.5th percentiles for each cytokine conc. 
    IL2_plus_top = np.percentile(IL2_plus, 97.5, axis=1)
    IL2_plus_bot = np.percentile(IL2_plus, 2.5, axis=1)
    IL2_minus_top = np.percentile(IL2_minus, 97.5, axis=1)
    IL2_minus_bot = np.percentile(IL2_minus, 2.5, axis=1)
    IL15_plus_top = np.percentile(IL15_plus, 97.5, axis=1)
    IL15_plus_bot = np.percentile(IL15_plus, 2.5, axis=1)
    IL15_minus_top = np.percentile(IL15_minus, 97.5, axis=1)
    IL15_minus_bot = np.percentile(IL15_minus, 2.5, axis=1)
    
    # plot shaded regions representing 95% confidence interval
    ax[0].fill_between(np.log10(cytokC), IL2_minus_top, IL2_minus_bot, color='darkorchid', alpha=0.5)
    ax[0].fill_between(np.log10(cytokC), IL15_minus_top, IL15_minus_bot, color='goldenrod', alpha=0.5)
    ax[1].fill_between(np.log10(cytokC), IL2_plus_top, IL2_plus_bot, color='darkorchid', alpha=0.5)
    ax[1].fill_between(np.log10(cytokC), IL15_plus_top, IL15_plus_bot, color='goldenrod', alpha=0.5)
       
    # plot experimental data
    path = os.path.dirname(os.path.abspath(__file__))
    data = pd.read_csv(join(path, "../data/IL2_IL15_extracted_data.csv")).values # imports file into pandas array
    ax[0].scatter(data[:,0], data[:,2], color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # IL2 in 2Ra-
    ax[0].scatter(data[:,0], data[:,3], color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # IL15 in 2Ra-
    ax[1].scatter(data[:,0], data[:,6], color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # IL2 in 2Ra+
    ax[1].scatter(data[:,0], data[:,7], color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # IL15 in 2Ra+
    ax[0].set(ylabel='Maximal p-STAT5 (% x 100)', xlabel='log10 of cytokine concentration (nM)', title='IL2Ra- YT-1 cells')
    ax[1].set(ylabel='Maximal p-STAT5 (% x 100)', xlabel='log10 of cytokine concentration (nM)', title='IL2Ra+ YT-1 cells')
    ax[0].legend()
    ax[1].legend()

    
def import_samples():
    """ This function imports the csv results into a numpy array called unkVec. """
    bmodel = build_model()
    n_params = nParams()

    path = os.path.dirname(os.path.abspath(__file__))
    trace = pm.backends.text.load(join(path, '../../IL2_model_results'), bmodel.M)
    kfwd = trace.get_values('kfwd', chains=[0])
    rxn = trace.get_values('rxn', chains=[0])
    endo_activeEndo = trace.get_values('endo', chains=[0])
    sortF = trace.get_values('sortF', chains=[0])
    kRec_kDeg = trace.get_values('kRec_kDeg', chains=[0])
    exprRates = trace.get_values('IL2Raexpr', chains=[0])
    
    unkVec = np.zeros((n_params, 500))
    for ii in range (0, 500):
        unkVec[:, ii] = np.array([0., 0., 0., 0., 0., 0., kfwd[ii], rxn[ii, 0], rxn[ii, 1], rxn[ii, 2], rxn[ii, 3], rxn[ii, 4], rxn[ii, 5], 1., 1., 1., 1., endo_activeEndo[ii, 0], endo_activeEndo[ii, 1], sortF[ii], kRec_kDeg[ii, 0], kRec_kDeg[ii, 1], exprRates[ii, 0], exprRates[ii, 1], exprRates[ii, 2], exprRates[ii, 3], 0., 0., 0., 0.])
    
    return unkVec

def violinPlots(ax, unkVec):
    """ Create violin plots of model posterior. """
    unkVec = unkVec.transpose()
    
    traf = pd.DataFrame(unkVec[:, 17:22])
    Rexpr = pd.DataFrame(unkVec[:, 22:26])
    
    traf.columns = traf_names()
    b = sns.violinplot(data=np.log10(traf), ax=ax[0])
    b.set_xticklabels(b.get_xticklabels(), rotation=40, rotation_mode="anchor", ha="right", fontsize=8, position=(0, 0.075))
    b.set_title("Trafficking parameters")
    b.set_ylabel("log10 of value")
    
    Rexpr.columns = ['IL2Ra', 'IL2Rb', 'gc', 'IL15Ra']
    c = sns.violinplot(data=np.log10(Rexpr), ax=ax[1])
    c.set_xticklabels(c.get_xticklabels(), rotation=40, rotation_mode="anchor", ha="right", fontsize=8, position=(0, 0.075))
    c.set_title("Receptor expression rates")
    c.set_ylabel("log10 of value")


def rateComp(ax, unkVec):
    """ This function compares the analogous reverse rxn distributions from IL2 and IL15 in a violin plot. """
    
    # assign values from unkVec
    kfwd, k4rev, k5rev, k16rev, k17rev, k22rev, k23rev = unkVec[6, :], unkVec[7, :], unkVec[8, :], unkVec[9, :], unkVec[10, :], unkVec[11, :], unkVec[12, :]
    
    # plug in values from measured constants
    kfbnd = 0.60 # Assuming on rate of 10^7 M-1 sec-1
    k1rev = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    k2rev = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    k13rev = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
    k14rev = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
    
    # make these scalar values arrays of size 500
    k1rev = np.full((500), k1rev)
    k2rev = np.full((500), k2rev)
    k13rev = np.full((500), k13rev)
    k14rev = np.full((500), k14rev)
    
    # detailed balance
    k10rev = 12.0 * k5rev / 1.5 # doi:10.1016/j.jmb.2004.04.038
    k11rev = 63.0 * k5rev / 1.5 # doi:10.1016/j.jmb.2004.04.038
    k12rev = k1rev * k11rev / k2rev # loop for IL2_IL2Ra_IL2Rb
    k9rev = k10rev * k11rev / k4rev
    k8rev = k10rev * k12rev / k5rev
    k24rev = k13rev * k23rev / k14rev
    k21rev = k22rev * k23rev / k16rev
    k20rev = k22rev * k24rev / k17rev
    
    # append analogous rxnrates for IL2/15
    k1_k13 = np.append(k1rev, k13rev)
    k2_k14 = np.append(k2rev, k14rev)
    k4_k16 = np.append(k4rev, k16rev)
    k5_k17 = np.append(k5rev, k17rev)
    k8_k20 = np.append(k8rev, k20rev)
    k9_k21 = np.append(k9rev, k21rev)
    k10_k22 = np.append(k10rev, k22rev)
    k11_k23 = np.append(k11rev, k23rev)
    k12_k24 = np.append(k12rev, k24rev)    
    
    # add each rate duo as separate column in dataframe
    df = pd.DataFrame({'k1_k13': k1_k13, 'k2_k14': k2_k14, 'k4_k16': k4_k16, 'k5_k17': k5_k17, 'k8_k20': k8_k20, 'k9_k21': k9_k21, 'k10_k22': k10_k22, 'k11_k23': k11_k23, 'k12_k24': k12_k24})
    
    # add labels for IL2 and IL15
    df['cytokine'] = 'IL2'
    df.loc[500:1000, 'cytokine'] = 'IL15'
    
    # melt into long form and take log value
    melted = pd.melt(df, id_vars='cytokine', var_name='rate', value_name='log10 of value')
    melted.loc[:, 'log10 of value'] = np.log10(melted.loc[:, 'log10 of value'])

    # plot with hue being cytokine species
    a = sns.violinplot(x='rate', y='log10 of value', data=melted, hue='cytokine', ax=ax)
    a.set_title("Analogous reverse reaction rates")
    
    