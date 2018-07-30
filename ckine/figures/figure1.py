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
    ax, f = getSetup((7, 6), (4, 4))

    subplotLabel(ax[0], 'A')
    unkVec = import_samples()
    #surf_perc(ax[0:4], 'IL2Rb', unkVec)
    #pstat_act(ax[4:6], unkVec)
    #surf_perc(ax[8:12], 'gc', unkVec)
    #violinPlots(ax[12:15], unkVec)
    rateComp(ax[15], unkVec)

    f.tight_layout()

    return f

def plot_structure(IL2vec, IL15vec, title, ax, x_axis, data_type, species='NONE'):
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
    # ax.legend()

def surf_perc(ax, species, unkVec):
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

    for ii in range(0,500):
        output = surf.calc(unkVec[:, ii], ts) * y_max
        IL2_1_plus = output[0:(size)]
        IL2_500_plus = output[(size):(size*2)]
        IL2_1_minus = output[(size*2):(size*3)]
        IL2_500_minus = output[(size*3):(size*4)]
        IL15_1_plus = output[(size*4):(size*5)]
        IL15_500_plus = output[(size*5):(size*6)]
        IL15_1_minus = output[(size*6):(size*7)]
        IL15_500_minus = output[(size*7):(size*8)]

        plot_structure(IL2_1_minus, IL15_1_minus, '1 nM and IL2Ra-', ax[0], ts, 'surf', species)
        plot_structure(IL2_500_minus, IL15_500_minus, "500 nM and IL2Ra-", ax[1], ts, 'surf', species)
        plot_structure(IL2_1_plus, IL15_1_plus, "1 nM and IL2Ra+", ax[2], ts, 'surf', species)
        plot_structure(IL2_500_plus, IL15_500_plus, "500 nM and IL2Ra+", ax[3], ts, 'surf', species)
        
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
    pstat5 = pstat()
    PTS = 30
    cytokC = np.logspace(-3.3, 2.7, PTS)
    y_max = 100.

    for ii in range(0,500):
        output = pstat5.calc(unkVec[:, ii], cytokC) * y_max
        IL2_plus = output[0:PTS]
        IL2_minus = output[PTS:(PTS*2)]
        IL15_plus = output[(PTS*2):(PTS*3)]
        IL15_minus = output[(PTS*3):(PTS*4)]

        plot_structure(IL2_minus, IL15_minus, "IL2Ra- YT-1 cells", ax[0], cytokC, 'act')
        plot_structure(IL2_plus, IL15_plus, "IL2Ra+ YT-1 cells", ax[1], cytokC, 'act')
       
    
    path = os.path.dirname(os.path.abspath(__file__))
    data = pd.read_csv(join(path, "../data/IL2_IL15_extracted_data.csv")).values # imports file into pandas array
    ax[0].scatter(data[:,0], data[:,2], color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # IL2 in 2Ra-
    ax[0].scatter(data[:,0], data[:,3], color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # IL15 in 2Ra-
    ax[1].scatter(data[:,0], data[:,6], color='darkorchid', marker='^', edgecolors='k', zorder=100, label='IL2') # IL2 in 2Ra+
    ax[1].scatter(data[:,0], data[:,7], color='goldenrod', marker='^', edgecolors='k', zorder=101, label='IL15') # IL15 in 2Ra+
    #ax[0].legend()
    #ax[1].legend()

    
def import_samples():
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
    
    rev_rxn = pd.DataFrame(unkVec[:, 7:13])
    traf = pd.DataFrame(unkVec[:, 17:22])
    Rexpr = pd.DataFrame(unkVec[:, 22:26])
    
    rev_rxn.columns = ['k4rev', 'k5rev', 'k16rev', 'k17rev', 'k22rev', 'k23rev']
    a = sns.violinplot(data=np.log10(rev_rxn), ax=ax[0])  # creates names based on dataframe columns
    a.set_xticklabels(a.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    
    traf.columns = traf_names()
    b = sns.violinplot(data=np.log10(traf), ax=ax[1])
    b.set_xticklabels(b.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))
    
    Rexpr.columns = ['IL2Ra', 'IL2Rb', 'gc', 'IL15Ra']
    c = sns.violinplot(data=np.log10(Rexpr), ax=ax[2])
    c.set_xticklabels(c.get_xticklabels(),
                       rotation=40,
                       rotation_mode="anchor",
                       ha="right",
                       fontsize=8,
                       position=(0, 0.075))


def rateComp(ax, unkVec):
    """ This function compares the analogous reverse rxn distributions from IL2 and IL15. """
    
    # assign values from unkVec
    kfwd, k4rev, k5rev, k16rev, k17rev, k22rev, k23rev = unkVec[6, :], unkVec[7, :], unkVec[8, :], unkVec[9, :], unkVec[10, :], unkVec[11, :], unkVec[12, :]
    
    # plug in values from measured constants
    kfbnd = 0.60 # Assuming on rate of 10^7 M-1 sec-1
    k1rev = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    k2rev = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    k13rev = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
    k14rev = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
    
    # detailed balance
    k10rev = 12.0 * k5rev / 1.5 # doi:10.1016/j.jmb.2004.04.038
    k11rev = 63.0 * k5rev / 1.5 # doi:10.1016/j.jmb.2004.04.038
    k12rev = k1rev * k11rev / k2rev # loop for IL2_IL2Ra_IL2Rb
    k9rev = k10rev * k11rev / k4rev
    k8rev = k10rev * k12rev / k5rev
    k24rev = k13rev * k23rev / k14rev
    k21rev = k22rev * k23rev / k16rev
    k20rev = k22rev * k24rev / k17rev
    
    k4_k16 = np.append(k4rev, k16rev)
    k5_k17 = np.append(k5rev, k17rev)
    k8_k20 = np.append(k8rev, k20rev)
    k9_k21 = np.append(k9rev, k21rev)
    k10_k22 = np.append(k10rev, k22rev)
    k11_k23 = np.append(k11rev, k23rev)
    k12_k24 = np.append(k12rev, k24rev)
    
    
    
    
    df = pd.DataFrame({'k4_k16': k4_k16, 'k5_k17': k5_k17, 'k8_k20': k8_k20, 'k9_k21': k9_k21, 'k10_k22': k10_k22, 'k11_k23': k11_k23, 'k12_k24': k12_k24})
    df['cytokine'] = 'IL2'
    df.loc[500:1000, 'cytokine'] = 'IL15'
    melted = pd.melt(df, id_vars='cytokine', var_name='rxn', value_name='vals')
    print(melted)
    melted.loc[:, 'vals'] = np.log10(melted.loc[:, 'vals'])
    print(melted)
    
    # test = pd.DataFrame({'rxn': ['k4/k16', 'k5/k17', 'k8/k20', 'k9/k21', 'k10/k22', 'k11/k23', 'k12/k24'], 'IL2': [k4rev, k5rev, k8rev, k9rev, k10rev, k11rev, k12rev], 'IL15': [k16rev, k17rev, k20rev, k21rev, k22rev, k23rev, k24rev] })
    # print(test)
    # test_melt = pd.melt(test, id_vars='rxn', var_name='cols',  value_name='vals')
    # print(test_melt)
    
    # A is cytokine type
    # B, C, D, are all the rxn rates
    
    # k1 vs k13
    # k2 vs k14
    # k4 vs k16
    # k5 vs k17
    # k8 vs k20
    # k9 vs k21
    # k10 vs k22
    # k11 vs k23
    # k12 vs k24
    
    # use pd.melt
    # id_vars will be IL2 and IL15
    # value_vars will be all the k* rxn rates
    a = sns.violinplot(x='rxn', y='vals', data=melted, hue='cytokine', ax=ax)
    
    