"""
This file includes various methods for flow cytometry analysis.
"""
import numpy as np
import pandas as pd
import seaborn as sns
import tensorly as tl
from tensorly.decomposition import non_negative_parafac, non_negative_parafac_hals
from copy import copy
from .figures.figureCommon import get_cellTypeDict

cellDict = get_cellTypeDict()


def makeTensor(sigDF, Variance=False):
    """Makes tensor of data with dimensions mutein x time point x concentration x cell type"""
    ligands = sigDF.Ligand.unique()
    tps = sigDF.Time.unique()
    concs = sigDF.Dose.unique()
    cellTypes = [
        "Treg $IL2Ra^{hi}$",
        "Treg",
        "Treg $IL2Ra^{lo}$",
        "Thelper $IL2Ra^{hi}$",
        "Thelper",
        "Thelper $IL2Ra^{lo}$",
        "CD8",
        "NK",
    ]
    tensor = np.empty((len(ligands), len(tps), len(concs), len(cellTypes)))
    tensor[:] = np.nan
    for i, lig in enumerate(ligands):
        for j, tp in enumerate(tps):
            for k, conc in enumerate(concs):
                for ii, cell in enumerate(cellTypes):
                    if not Variance:
                        entry = sigDF.loc[
                            (sigDF.Ligand == lig)
                            & (sigDF.Time == tp)
                            & (sigDF.Dose == conc)
                            & (sigDF.Cell == cell)
                        ].Mean.values
                    else:
                        entry = sigDF.loc[
                            (sigDF.Ligand == lig)
                            & (sigDF.Time == tp)
                            & (sigDF.Dose == conc)
                            & (sigDF.Cell == cell)
                        ].Variance.values
                    if len(entry) >= 1:
                        tensor[i, j, k, ii] = np.mean(entry)
    # Normalize
    for i, _ in enumerate(cellTypes):
        tensor[:, :, :, i][~np.isnan(tensor[:, :, :, i])] /= np.nanmax(
            tensor[:, :, :, i]
        )

    return tensor


def factorTensor(tensor, numComps):
    """Takes Tensor, and mask and returns tensor factorized form."""
    tfac = non_negative_parafac(
        np.nan_to_num(tensor),
        rank=numComps,
        mask=np.isfinite(tensor),
        n_iter_max=5000,
        tol=1e-9,
    )
    tensor = tensor.copy()
    tensor[np.isnan(tensor)] = tl.cp_to_tensor(tfac)[np.isnan(tensor)]
    return non_negative_parafac_hals(tensor, numComps, n_iter_max=5000)


def R2Xplot(ax, tensor, compNum):
    """Creates R2X plot for non-neg CP tensor decomposition"""
    varHold = np.zeros(compNum)
    for i in range(1, compNum + 1):
        tFac = factorTensor(tensor, i)
        varHold[i - 1] = calcR2X(tensor, tFac)

    ax.scatter(np.arange(1, compNum + 1), varHold, c="k", s=20.0)
    ax.set(
        title="R2X",
        ylabel="Variance Explained",
        xlabel="Number of Components",
        ylim=(0, 1),
        xlim=(0, compNum + 0.5),
        xticks=np.arange(0, compNum + 1),
    )


def calcR2X(tensorIn, tensorFac):
    """Calculate R2X."""
    tErr = np.nanvar(tl.cp_to_tensor(tensorFac) - tensorIn)
    return 1.0 - tErr / np.nanvar(tensorIn)


def plot_tFac_Ligs(ax, tFac, respDF, numComps=3):
    """Plots tensor factorization of Ligands"""
    ligands = respDF.Ligand.unique()
    mutFacs = tFac[1][0]
    tFacDFlig = pd.DataFrame()
    for i in range(0, numComps):
        tFacDFlig = pd.concat(
            [
                tFacDFlig,
                pd.DataFrame(
                    {
                        "Component_Val": mutFacs[:, i],
                        "Component": (i + 1),
                        "Ligand": ligands,
                    }
                ),
            ]
        )

    tFacDFlig["Valency"], tFacDFlig["Lig Name"] = 0, 0
    tFacDFlig = tFacDFlig.reset_index()
    for index, row in tFacDFlig.iterrows():
        if row.Ligand.split(" ")[-1] == "(Mono)":
            valency = 1
        elif row.Ligand.split(" ")[-1] == "(Biv)":
            valency = 2
        tFacDFlig.iloc[index, tFacDFlig.columns.get_loc("Valency")] = valency
        tFacDFlig.iloc[index, tFacDFlig.columns.get_loc("Lig Name")] = " ".join(
            row.Ligand.split(" ")[0:-1]
        )

    sns.barplot(data=tFacDFlig, ax=ax, x="Component_Val", y="Ligand", hue="Component")
    ax.set(xlabel="Component Weight")
    return tFacDFlig


def plot_tFac_Time(ax, tFac, respDF):
    """Plots tensor factorization of tps"""
    tps = respDF.Time.unique()
    timeFacs = tFac[1][1]
    tps = np.append(tps, [0])
    timeFacs = np.append(timeFacs, np.zeros((1, timeFacs.shape[1])), axis=0)

    markersTimes = [".", ".", ".", "."]
    for i in range(0, timeFacs.shape[1]):
        ax.plot(
            tps,
            timeFacs[:, i],
            marker=markersTimes[i],
            label="Component " + str(i + 1),
            markersize=5,
        )

    ax.set(
        title="Time",
        xlabel="Time (hrs)",
        xlim=(0, 4),
        ylabel="Component Weight",
        ylim=(0, 1),
    )


def plot_tFac_Conc(ax, tFac, respDF):
    """Plots tensor factorization of Conc"""
    concs = respDF.Dose.unique()
    concFacs = tFac[1][2]

    markersConcs = [".", ".", ".", "."]
    for i in range(0, concFacs.shape[1]):
        ax.plot(
            concs,
            concFacs[:, i],
            marker=markersConcs[i],
            label="Component " + str(i + 1),
            markersize=5,
        )

    ax.set(
        title="Concentration",
        xlabel="Concentration (nM)",
        xlim=(concs[-1], concs[0]),
        ylabel="Component Weight",
        ylim=(0, 1),
        xscale="log",
    )


def plot_tFac_Cells(ax, tFac, respDF, numComps=3):
    """Plots tensor factorization of cells"""
    cells = [
        "Treg $IL2Ra^{hi}$",
        "Treg",
        "Treg $IL2Ra^{lo}$",
        "Thelper $IL2Ra^{hi}$",
        "Thelper",
        "Thelper $IL2Ra^{lo}$",
        "CD8",
        "NK",
    ]
    cellFacs = tFac[1][3]
    tFacDFcell = pd.DataFrame()
    for i in range(0, numComps):
        tFacDFcell = pd.concat(
            [
                tFacDFcell,
                pd.DataFrame(
                    {
                        "Component_Val": cellFacs[:, i],
                        "Component": (i + 1),
                        "Cell": cells,
                    }
                ),
            ]
        )

    tFacDFcell = tFacDFcell.replace(cellDict)
    sns.barplot(data=tFacDFcell, ax=ax, x="Cell", y="Component_Val", hue="Component")
    ax.set(ylabel="Component Weight")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


def facScatterPlot(ax, tFacDFLig):
    """Plots tensor factorization of cells"""
    scattDF = copy(tFacDFLig.drop("index", axis=1))
    for ligand in scattDF.Ligand.unique():
        for valency in scattDF.loc[scattDF.Ligand == ligand].Valency.unique():
            isoDF = scattDF.loc[
                (scattDF.Ligand == ligand) & (scattDF.Valency == valency)
            ]
            scattDF = pd.concat(
                [
                    scattDF,
                    pd.DataFrame(
                        {
                            "Component 1 + 3": isoDF.loc[
                                isoDF.Component == 1
                            ].Component_Val.values
                            + isoDF.loc[isoDF.Component == 3].Component_Val.values,
                            "Component 2": isoDF.loc[
                                isoDF.Component == 2
                            ].Component_Val.values,
                            "Valency": valency,
                        }
                    ),
                ]
            )
    sns.scatterplot(
        data=scattDF,
        x="Component 1 + 3",
        y="Component 2",
        hue="Valency",
        style="Valency",
        size="Valency",
        ax=ax,
    )
    ax.set(xlim=(0, 0.8), ylim=(0, 0.8))
