"""
This file includes various methods for flow cytometry analysis.
"""
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import tensorly as tl
import matplotlib.cm as cm
from tensorly.decomposition import non_negative_parafac


def makeTensor(sigDF, Variance=False):
    """Makes tensor of data with dimensions mutein x time point x concentration x cell type"""
    ligands = sigDF.Ligand.unique()
    tps = sigDF.Time.unique()
    concs = sigDF.Dose.unique()
    cellTypes = sigDF.Cell.unique()
    tensor = np.empty((len(ligands), len(tps), len(concs), len(cellTypes)))
    tensor[:] = np.nan
    for i, lig in enumerate(ligands):
        for j, tp in enumerate(tps):
            for k, conc in enumerate(concs):
                for ii, cell in enumerate(cellTypes):
                    if not Variance:
                        entry = sigDF.loc[(sigDF.Ligand == lig) & (sigDF.Time == tp) & (sigDF.Dose == conc) & (sigDF.Cell == cell)].Mean.values
                    else:
                        entry = sigDF.loc[(sigDF.Ligand == lig) & (sigDF.Time == tp) & (sigDF.Dose == conc) & (sigDF.Cell == cell)].Variance.values
                    if len(entry) >= 1:
                        tensor[i, j, k, ii] = np.mean(entry)
    # Normalize
    for i, _ in enumerate(cellTypes):
        tensor[:, :, :, i][~np.isnan(tensor[:, :, :, i])] /= np.nanmax(tensor[:, :, :, i])

    return tensor


def factorTensor(tensor, numComps):
    """ Takes Tensor, and mask and returns tensor factorized form. """
    return non_negative_parafac(np.nan_to_num(tensor), rank=numComps, mask=np.isfinite(tensor), n_iter_max=5000, init='svd', random_state=0)


def R2Xplot(ax, tensor, compNum):
    """Creates R2X plot for non-neg CP tensor decomposition"""
    varHold = np.zeros(compNum)
    for i in range(1, compNum + 1):
        tFac = factorTensor(tensor, i)
        varHold[i - 1] = calcR2X(tensor, tFac)

    ax.scatter(np.arange(1, compNum + 1), varHold, c='k', s=20.)
    ax.set(title="R2X", ylabel="Variance Explained", xlabel="Number of Components", ylim=(0, 1), xlim=(0, compNum + 1), xticks=np.arange(0, compNum + 1))


def calcR2X(tensorIn, tensorFac):
    """ Calculate R2X. """
    tErr = np.nanvar(tl.cp_to_tensor(tensorFac) - tensorIn)
    return 1.0 - tErr / np.nanvar(tensorIn)


def plot_tFac_Ligs(ax, tFac, respDF):
    """Plots tensor factorization of Ligands"""
    ligands = respDF.Ligand.unique()
    mutFacs = tFac[1][0]
    tFacDFlig = pd.DataFrame({"Cmpnt. 1": mutFacs[:, 0], "Cmpnt. 2": mutFacs[:, 1], "Cmpnt. 3": mutFacs[:, 2]}, index=ligands)

    ax.tick_params(axis='x', labelrotation=45)
    sns.heatmap(tFacDFlig, ax=ax, vmin=0., vmax=1.0, cbar_kws={'label': 'Component Value'})


def plot_tFac_Time(ax, tFac, respDF):
    """Plots tensor factorization of tps"""
    tps = respDF.Time.unique()
    timeFacs = tFac[1][1]

    markersTimes = ["^", "*", "D"]
    for i in range(0, timeFacs.shape[1]):
        ax.plot(tps, timeFacs[:, i], marker=markersTimes[i], label="Component " + str(i + 1))

    ax.legend()
    ax.set(title="Time", xlabel="Time (hrs)", xlim=(0, 4), ylabel="Component", ylim=(0, 1))


def plot_tFac_Conc(ax, tFac, respDF):
    """Plots tensor factorization of Conc"""
    concs = respDF.Dose.unique()
    concFacs = tFac[1][2]

    markersConcs = ["^", "*", "D"]
    for i in range(0, concFacs.shape[1]):
        ax.plot(concs, concFacs[:, i], marker=markersConcs[i], label="Component " + str(i + 1))

    ax.legend()
    ax.set(title="Concentration", xlabel="Concentration (nM)", xlim=(concs[-1], concs[0]), ylabel="Component", ylim=(0, 1), xscale='log')


def plot_tFac_Cells(ax, tFac, respDF):
    """Plots tensor factorization of cells"""
    cells = respDF.Cell.unique()
    cellFacs = tFac[1][3]
    tFacDFcell = pd.DataFrame({"Cmpnt. 1": cellFacs[:, 0], "Cmpnt. 2": cellFacs[:, 1], "Cmpnt. 3": cellFacs[:, 2]}, index=cells)
    ax.tick_params(axis='x', labelrotation=45)

    sns.heatmap(tFacDFcell, ax=ax, vmin=0., vmax=1.0, cbar_kws={'label': 'Component Value'})
