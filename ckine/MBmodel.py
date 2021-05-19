"""
Implementation of a simple multivalent binding model.
"""

from os.path import dirname, join
import numpy as np
import pandas as pd
from scipy.optimize import root
from .imports import import_pstat_all, getBindDict, importReceptors

path_here = dirname(dirname(__file__))
KxStarP = 1.126e-12


def getKxStar():
    return KxStarP


def Req_func2(Req, L0, KxStar, Rtot, Cplx, Ctheta, Kav):
    Psi = Req * Kav * KxStar
    Psi = np.pad(Psi, ((0, 0), (0, 1)), constant_values=1)
    Psirs = np.sum(Psi, axis=1).reshape(-1, 1)
    Psinorm = (Psi / Psirs)[:, :-1]

    Rbound = L0 / KxStar * np.sum(Ctheta.reshape(-1, 1) * np.dot(Cplx, Psinorm) * np.exp(np.dot(Cplx, np.log1p(Psirs - 1))), axis=0)
    return Req + Rbound - Rtot


def polyc(L0, KxStar, Rtot, Cplx, Ctheta, Kav):
    """
    The main function to be called for multivalent binding
    :param L0: concentration of ligand complexes
    :param KxStar: Kx for detailed balance correction
    :param Rtot: numbers of each receptor on the cell
    :param Cplx: the monomer ligand composition of each complex
    :param Ctheta: the composition of complexes
    :param Kav: Ka for monomer ligand to receptors
    :return:
        Lbound: a list of Lbound of each complex
        Rbound: a list of Rbound of each kind of receptor
    """
    # Consistency check
    Kav = np.array(Kav)
    assert Kav.ndim == 2
    Rtot = np.array(Rtot, dtype=np.float)
    assert Rtot.ndim == 1
    Cplx = np.array(Cplx)
    assert Cplx.ndim == 2
    Ctheta = np.array(Ctheta)
    assert Ctheta.ndim == 1

    assert Kav.shape[0] == Cplx.shape[1]
    assert Kav.shape[1] == Rtot.size
    assert Cplx.shape[0] == Ctheta.size
    Ctheta = Ctheta / np.sum(Ctheta)

    # Solve Req
    lsq = root(Req_func2, Rtot, method="lm", args=(L0, KxStar, Rtot, Cplx, Ctheta, Kav), options={"maxiter": 3000})
    assert lsq["success"], "Failure in rootfinding. " + str(lsq)
    Req = lsq["x"].reshape(1, -1)

    # Calculate the results
    Psi = np.ones((Kav.shape[0], Kav.shape[1] + 1))
    Psi[:, : Kav.shape[1]] *= Req * Kav * KxStar
    Psirs = np.sum(Psi, axis=1).reshape(-1, 1)
    Psinorm = (Psi / Psirs)[:, :-1]

    Rbound = L0 / KxStar * Ctheta.reshape(-1, 1) * np.dot(Cplx, Psinorm) * np.exp(np.dot(Cplx, np.log1p(Psirs - 1)))
    assert Rbound.shape[0] == len(Ctheta)
    assert Rbound.shape[1] == len(Rtot)
    return Rbound


def cytBindingModel(mut, val, doseVec, cellType, x=False, date=False, binNum=False):
    """Runs binding model for a given mutein, valency, dose, and cell type."""
    recQuantDF = importReceptors()
    recCount = np.ravel([recQuantDF.loc[(recQuantDF.Receptor == "IL2Ra") & (recQuantDF["Cell Type"] == cellType)].Mean.values,
                             recQuantDF.loc[(recQuantDF.Receptor == "IL2Rb") & (recQuantDF["Cell Type"] == cellType)].Mean.values])

    mutAffDF = pd.read_csv(join(path_here, "ckine/data/WTmutAffData.csv"))
    Affs = mutAffDF.loc[(mutAffDF.Mutein == mut)]
    Affs = np.power(np.array([Affs["IL2RaKD"].values, Affs["IL2RBGKD"].values]) / 1e9, -1)
    Affs = np.reshape(Affs, (1, -1))
    Affs = np.repeat(Affs, 2, axis=0)
    np.fill_diagonal(Affs, 1e2)  # Each cytokine can only bind one a and one b

    if doseVec.size == 1:
        doseVec = np.array([doseVec])
    output = np.zeros(doseVec.size)

    for i, dose in enumerate(doseVec):
        if x:
            output[i] = polyc(dose / 1e9, np.power(10, x[0]), recCount, [[val, val]], [1.0], Affs)[0][1]
        else:
            output[i] = polyc(dose / 1e9, KxStarP, recCount, [[val, val]], [1.0], Affs)[0][1]  # IL2RB binding only
    if date:
        convDict = getBindDict()
        if cellType[-1] == "$": #if it is a binned pop, use ave fit
            output *= convDict.loc[(convDict.Date == date) & (convDict.Cell == cellType[0:-13])].Scale.values
        else:
            output *= convDict.loc[(convDict.Date == date) & (convDict.Cell == cellType)].Scale.values
    return output


def runFullModel(x=False, time=[0.5], saveDict=False, singleCell=False):
    """Runs model for all data points and outputs date conversion dict for binding to pSTAT. Can be used to fit Kx"""
    statDF = import_pstat_all(singleCell)
    statDF = statDF.loc[(statDF.Ligand != "H16L N-term (Mono)") & (statDF.Ligand != "IL15 (Mono)")]
    statDF = statDF.loc[(statDF.Time.isin(time))]

    if x:
        # If Minimizing, scale all so that cells are weighted identically
        statDF.loc[(statDF.Cell == "Treg"), "Mean"] *= 1000 / statDF.loc[(statDF.Cell == "Treg"), "Mean"].mean()
        statDF.loc[(statDF.Cell == "Thelper"), "Mean"] *= 1000 / statDF.loc[(statDF.Cell == "Thelper"), "Mean"].mean()
        statDF.loc[(statDF.Cell == "CD8"), "Mean"] *= 1000 / statDF.loc[(statDF.Cell == "CD8"), "Mean"].mean()
        statDF.loc[(statDF.Cell == "NK"), "Mean"] *= 1000 / statDF.loc[(statDF.Cell == "NK"), "Mean"].mean()

    dateConvDF = pd.DataFrame(columns={"Date", "Scale", "Cell"})
    masterSTAT = pd.DataFrame(columns={"Ligand", "Date", "Cell", "Time", "Dose", "Valency", "Experimental", "Predicted"})
    dates = statDF.Date.unique()

    for (date, lig, conc, cell, time), group in statDF.groupby(["Date", "Ligand", "Dose", "Cell", "Time"]):
        if lig[-5::] == "(Biv)":
            val = 2
            ligName = lig[0:-6]
        else:
            val = 1
            ligName = lig[0:-7]

            entry = group.Mean.values
            if len(entry) >= 1:
                expVal = np.mean(entry)
                predVal = cytBindingModel(ligName, val, conc, cell, x)
                masterSTAT = masterSTAT.append(pd.DataFrame({"Ligand": ligName, "Date": date, "Cell": cell, "Dose": conc,
                                                             "Time": time, "Valency": val, "Experimental": expVal, "Predicted": predVal}))
        for date in dates:
            for cell in masterSTAT.Cell.unique():
                if cell[-1] == "$": #if it is a binned pop, use ave fit
                    predVecBin = masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell)].Predicted.values

                    expVec = masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell[0:-13])].Experimental.values
                    predVec = masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell[0:-13])].Predicted.values
                    slope = np.linalg.lstsq(np.reshape(predVec, (-1, 1)), np.reshape(expVec, (-1, 1)), rcond=None)[0][0]

                    masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell), "Predicted"] = predVecBin * slope
                else:
                    expVec = masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell)].Experimental.values
                    predVec = masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell)].Predicted.values
                    slope = np.linalg.lstsq(np.reshape(predVec, (-1, 1)), np.reshape(expVec, (-1, 1)), rcond=None)[0][0]
                    masterSTAT.loc[(masterSTAT.Date == date) & (masterSTAT.Cell == cell), "Predicted"] = predVec * slope
                    dateConvDF = dateConvDF.append(pd.DataFrame({"Date": date, "Scale": slope, "Cell": cell}))
    if saveDict:
        dateConvDF.set_index("Date").to_csv(join(path_here, "ckine/data/BindingConvDict.csv"))

    if x:
        return np.linalg.norm(masterSTAT.Predicted.values - masterSTAT.Experimental.values)
    else:
        return masterSTAT
