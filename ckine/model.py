"""
A file that includes the model and important helper functions.
"""
import os
import ctypes as ct
import numpy as np
from collections import OrderedDict


filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "./ckine.so")
libb = ct.cdll.LoadLibrary(filename)
pcd = ct.POINTER(ct.c_double)
libb.fullModel_C.argtypes = (pcd, ct.c_double, pcd, pcd)
libb.runCkine.argtypes = (pcd, ct.c_uint, pcd, pcd, ct.c_bool, ct.c_double, pcd)
libb.runCkineParallel.argtypes = (pcd, pcd, ct.c_uint, ct.c_uint, pcd, ct.c_double, pcd)
libb.runCkineSParallel.argtypes = (pcd, pcd, ct.c_uint, ct.c_uint, pcd, pcd, pcd, ct.c_double, pcd)

__nSpecies = 62


def nSpecies():
    """ Returns the total number of species in the model. """
    return __nSpecies


__halfL = 28


def halfL():
    """ Returns the number of species on the surface alone. """
    return __halfL


__nParams = 30
__rxParams = 60


def nParams():
    """ Returns the length of the rxntfR vector. """
    return __nParams


def rxParams():
    """ Returns the length of the rxntfR vector. """
    return __rxParams


__internalStrength = 0.5  # strength of endosomal activity relative to surface


def internalStrength():
    """Returns the internalStrength of endosomal activity."""
    return __internalStrength


__internalV = 623.0  # endosomal volume


def internalV():
    """ Returns __internalV. """
    return __internalV


__nRxn = 17


def nRxn():
    """ Returns the length of the rxn rates vector (doesn't include traf rates). """
    return __nRxn


def runCkineU_IL2(tps, rxntfr):
    """ Standard version of solver that returns species abundances given times and unknown rates. """
    rxntfr = rxntfr.copy()
    assert rxntfr.size == 15
    assert np.all(rxntfr >= 0.0)

    yOut = np.zeros((tps.size, __nSpecies), dtype=np.float64)

    rxntfr = getRateVec(rxntfr)

    retVal = libb.runCkine(tps.ctypes.data_as(ct.POINTER(ct.c_double)), tps.size, yOut.ctypes.data_as(ct.POINTER(ct.c_double)), rxntfr.ctypes.data_as(ct.POINTER(ct.c_double)), True, 0.0, None)

    assert retVal >= 0  # make sure solver worked

    return yOut


def runIL2simple(unkVec, input_params, IL, CD25=1.0, tps=None, input_receptors=None, ligandDegradation=False):
    """ Version to focus on IL2Ra/Rb affinity adjustment. """

    if tps is None:
        tps = np.array([500.0])

    kfwd, k4rev, k5rev = unkVec[6], unkVec[7], unkVec[8]

    k1rev = 0.6 * 10 * input_params[0]
    k2rev = 0.6 * 144 * input_params[1]
    k11rev = 63.0 * k5rev / 1.5 * input_params[1]

    if input_receptors is not None:
        IL2Ra = input_receptors[0] * CD25
        IL2Rb = input_receptors[1]
        gc = input_receptors[2]
    else:
        IL2Ra, IL2Rb, gc = unkVec[22] * CD25, unkVec[23], unkVec[24]

    # IL, kfwd, k1rev, k2rev, k4rev, k5rev, k11rev, R, R, R
    rxntfr = np.array([IL, kfwd, k1rev, k2rev, k4rev, k5rev, k11rev, IL2Ra, IL2Rb, gc, k1rev * input_params[2], k2rev * input_params[2],
                       k4rev * input_params[2], k5rev * input_params[2], k11rev * input_params[2]])  # input_params[2] represents endosomal binding affinity relative to surface affinity

    yOut = runCkineU_IL2(tps, rxntfr)

    if ligandDegradation:
        # rate of ligand degradation
        return ligandDeg(yOut[0], sortF=unkVec[19], kDeg=unkVec[21], cytokineIDX=0)

    active_ckine = np.zeros(yOut.shape[0])

    # calculate for each time point
    for i in range(yOut.shape[0]):
        active_ckine[i] = getTotalActiveCytokine(0, yOut[i, :])

    return active_ckine


def runCkineU(tps, rxntfr, preT=0.0, prestim=None):
    """ Standard version of solver that returns species abundances given times and unknown rates. """
    return runCkineUP(tps, np.atleast_2d(rxntfr.copy()), preT, prestim)


def runCkineUP(tps, rxntfr, preT=0.0, prestim=None):
    """ Version of runCkine that runs in parallel. """
    tps = np.array(tps)
    assert rxntfr.size % __nParams == 0
    assert rxntfr.shape[1] == __nParams

    assert (rxntfr[:, 19] < 1.0).all()  # Check that sortF won't throw
    assert np.all(np.any(rxntfr > 0.0, axis=1))  # make sure at least one element is >0 for all rows

    yOut = np.zeros((rxntfr.shape[0] * tps.size, __nSpecies), dtype=np.float64)

    rxntfr = getRateVec(rxntfr)

    if preT != 0.0:
        assert preT > 0.0
        assert prestim.size == 6
        prestim = prestim.ctypes.data_as(ct.POINTER(ct.c_double))

    retVal = libb.runCkineParallel(
        rxntfr.ctypes.data_as(ct.POINTER(ct.c_double)), tps.ctypes.data_as(ct.POINTER(ct.c_double)), tps.size, rxntfr.shape[0], yOut.ctypes.data_as(ct.POINTER(ct.c_double)), preT, prestim
    )

    assert retVal >= 0  # make sure solver worked

    return yOut


def runCkineSP(tps, rxntfr, actV, preT=0.0, prestim=None):
    """ Version of runCkine that runs in parallel. """
    tps = np.array(tps)
    assert rxntfr.size % __nParams == 0
    assert rxntfr.shape[1] == __nParams
    assert (rxntfr[:, 19] < 1.0).all()  # Check that sortF won't throw

    yOut = np.zeros((rxntfr.shape[0] * tps.size), dtype=np.float64)

    rxntfr = getRateVec(rxntfr)
    sensV = np.zeros((rxntfr.shape[0] * tps.size, __rxParams), dtype=np.float64, order="C")

    if preT != 0.0:
        assert preT > 0.0
        assert prestim.size == 6
        prestim = prestim.ctypes.data_as(ct.POINTER(ct.c_double))

    retVal = libb.runCkineSParallel(
        rxntfr.ctypes.data_as(ct.POINTER(ct.c_double)),
        tps.ctypes.data_as(ct.POINTER(ct.c_double)),
        tps.size,
        rxntfr.shape[0],
        yOut.ctypes.data_as(ct.POINTER(ct.c_double)),
        sensV.ctypes.data_as(ct.POINTER(ct.c_double)),
        actV.ctypes.data_as(ct.POINTER(ct.c_double)),
        preT,
        prestim,
    )

    sensV = condenseSENV(sensV)

    return (yOut, retVal, sensV)


def fullModel(y, t, rxntfr):
    """ Implement the full model based on dydt, trafficking, expression. """
    assert rxntfr.size == __nParams

    yOut = np.zeros_like(y)

    libb.fullModel_C(y.ctypes.data_as(ct.POINTER(ct.c_double)), t, yOut.ctypes.data_as(ct.POINTER(ct.c_double)), rxntfr.ctypes.data_as(ct.POINTER(ct.c_double)))

    return yOut


__active_species_IDX = np.zeros(__halfL, dtype=np.bool)
__active_species_IDX[np.array([7, 8, 14, 15, 18, 21, 24, 27])] = 1


def getActiveSpecies():
    """ Return a vector that indicates which species are active. """
    return __active_species_IDX


def getTotalActiveSpecies():
    """ Return a vector of all the species (surface + endosome) which are active. """
    activity = getActiveSpecies()
    return np.concatenate((activity, __internalStrength * activity, np.zeros(6)))


def getCytokineSpecies():
    """ Returns a list of vectors for which species are bound to which cytokines. """
    return list((np.arange(3, 9), np.arange(10, 16), np.arange(17, 19), np.arange(20, 22), np.arange(23, 25), np.arange(26, 28)))


def getSurfaceIL2RbSpecies():
    """ Returns a list of vectors for which surface species contain the IL2Rb receptor. """
    condense = np.zeros(__nSpecies)
    condense[np.array([1, 4, 5, 7, 8, 11, 12, 14, 15])] = 1
    return condense


def getSurfaceGCSpecies():
    """ Returns a list of vectors for which surface species contain the gc receptor. """
    condense = np.zeros(__nSpecies)
    condense[np.array([2, 6, 7, 8, 13, 14, 15, 18, 21])] = 1
    return condense


def getActiveCytokine(cytokineIDX, yVec):
    """ Get amount of active species. """
    return ((yVec * getActiveSpecies())[getCytokineSpecies()[cytokineIDX]]).sum()


def getTotalActiveCytokine(cytokineIDX, yVec):
    """ Get amount of surface and endosomal active species. """
    assert yVec.ndim == 1
    return getActiveCytokine(cytokineIDX, yVec[0:__halfL]) + __internalStrength * getActiveCytokine(cytokineIDX, yVec[__halfL: __halfL * 2])


def surfaceReceptors(y):
    """This function takes in a vector y and returns the amounts of the 8 surface receptors"""
    IL2Ra = np.sum(y[np.array([0, 3, 5, 6, 8])])
    IL2Rb = np.sum(y[np.array([1, 4, 5, 7, 8, 11, 12, 14, 15])])
    gc = np.sum(y[np.array([2, 6, 7, 8, 13, 14, 15, 18, 21])])
    IL15Ra = np.sum(y[np.array([9, 10, 12, 13, 15])])
    IL7Ra = np.sum(y[np.array([16, 17, 18])])
    IL9R = np.sum(y[np.array([19, 20, 21])])
    IL4Ra = np.sum(y[np.array([22, 23, 24])])
    IL21Ra = np.sum(y[np.array([25, 26, 27])])
    return np.array([IL2Ra, IL2Rb, gc, IL15Ra, IL7Ra, IL9R, IL4Ra, IL21Ra])


def totalReceptors(yVec):
    """This function takes in a vector y and returns the amounts of all 8 receptors in both cell compartments"""
    return surfaceReceptors(yVec) + __internalStrength * surfaceReceptors(yVec[__halfL: __halfL * 2])


def ligandDeg(yVec, sortF, kDeg, cytokineIDX):
    """ This function calculates rate of total ligand degradation. """
    yVec_endo_species = yVec[__halfL: (__halfL * 2)].copy()  # get all endosomal complexes
    yVec_endo_lig = yVec[(__halfL * 2)::].copy()  # get all endosomal ligands
    sum_active = np.sum(getActiveCytokine(cytokineIDX, yVec_endo_species))
    __cytok_species_IDX = np.zeros(__halfL, dtype=np.bool)  # create array of size halfL
    __cytok_species_IDX[getCytokineSpecies()[cytokineIDX]] = 1  # assign 1's for species corresponding to the cytokineIDX
    sum_total = np.sum(yVec_endo_species * __cytok_species_IDX)
    sum_inactive = (sum_total - sum_active) * sortF  # scale the inactive species by sortF
    return kDeg * (((sum_inactive + sum_active) * __internalStrength) + (yVec_endo_lig[cytokineIDX] * __internalV))  # can assume all free ligand and active species are degraded at rate kDeg


def receptor_expression(receptor_abundance, endo, kRec, sortF, kDeg):
    """ Uses receptor abundance (from flow) and trafficking rates to calculate receptor expression rate at steady state. """
    rec_ex = (receptor_abundance * endo) / (1.0 + ((kRec * (1.0 - sortF)) / (kDeg * sortF)))
    return rec_ex


def condenseSENV(sensVin):
    """ Condense sensitivities down into the old rxnRates format. """
    sensVin[:, 7:27] += sensVin[:, 27:47] * 5.0
    sensVin[:, 10] += 12.0 * sensVin[:, 11] / 1.5 + 63.0 * sensVin[:, 12] / 1.5
    sensV = sensVin[:, np.array([0, 1, 2, 3, 4, 5, 6, 9, 10, 15, 16, 17, 18, 20, 22, 24, 26, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59])]

    return sensV


def getparamsdict(rxntfr):
    """Where rate vectors and constants are defined, organized in an ordered dictionary"""
    rd = OrderedDict()
    modtype = rxntfr.size

    if modtype == 15:
        rd['IL2'] = rxntfr[0]
        rd['IL15'], rd['IL7'], rd['IL9'], rd['IL4'], rd['IL21'] = 0, 0, 0, 0, 0
        rd['kfwd'] = rxntfr[1]
        rd['surface.k1rev'] = rxntfr[2]
        rd['surface.k2rev'] = rxntfr[3]
        rd['surface.k4rev'] = rxntfr[4]
        rd['surface.k5rev'] = rxntfr[5]  # 9 #10
        rd['surface.k10rev'] = 12.0 * rd['surface.k5rev'] / 1.5
        rd['surface.k11rev'] = rxntfr[6]
        rd['surface.k13rev'] = 1.0
        rd['surface.k14rev'] = 1.0
        rd['surface.k16rev'] = 1.0
        rd['surface.k17rev'] = 1.0
        rd['surface.k22rev'] = 1.0
        rd['surface.k23rev'] = 1.0
        rd['surface.k25rev'] = 1.0
        rd['surface.k27rev'] = 1.0
        rd['surface.k29rev'] = 1.0
        rd['surface.k31rev'] = 1.0
        rd['surface.k32rev'] = 1.0
        rd['surface.k33rev'] = 1.0
        rd['surface.k34rev'] = 1.0
        rd['surface.k35rev'] = 1.0
        rd['endosome.k1rev'] = rxntfr[10]
        rd['endosome.k2rev'] = rxntfr[11]
        rd['endosome.k4rev'] = rxntfr[12]
        rd['endosome.k5rev'] = rxntfr[13]
        rd['endosome.k10rev'] = 12 * rd['endosome.k5rev'] / 1.5
        rd['endosome.k11rev'] = rxntfr[14]

        for ii in ('13', '14', '16', '17', '22', '23', '25', '27', '29', '31', '32', '33', '34', '35'):
            rd['endosome.k' + ii + 'rev'] = rd['surface.k' + ii + 'rev']

        rd['endo'] = 0.08221
        rd['activeEndo'] = 2.52654
        rd['sortF'] = 0.16024
        rd['kRec'] = 0.10017
        rd['kDeg'] = 0.00807

        rd['Rexpr_2Ra'], rd['Rexpr_2Rb'], rd['Rexpr_gc'], = rxntfr[7:10]
        rd['Rexpr_15Ra'], rd['Rexpr_7R'], rd['Rexpr_9R'], rd['Rexpr_4Ra'], rd['Rexpr_21Ra'] = 0, 0, 0, 0, 0

    else:
        kfbnd = 0.60
        rd['IL2'], rd['IL15'], rd['IL7'], rd['IL9'], rd['IL4'], rd['IL21'], rd['kfwd'] = tuple(rxntfr[0:7])
        rd['surface.k1rev'] = kfbnd * 10.0  # 7
        rd['surface.k2rev'] = kfbnd * 144.0
        rd['surface.k4rev'], rd['surface.k5rev'] = rxntfr[7], rxntfr[8]  # 9 #10
        rd['surface.k10rev'] = 12.0 * rd['surface.k5rev'] / 1.5
        rd['surface.k11rev'] = 63.0 * rd['surface.k5rev'] / 1.5
        rd['surface.k13rev'] = kfbnd * 0.065
        rd['surface.k14rev'] = kfbnd * 438.0
        rd['surface.k16rev'] = rxntfr[9]
        rd['surface.k17rev'] = rxntfr[10]  # 16
        rd['surface.k22rev'] = rxntfr[11]
        rd['surface.k23rev'] = rxntfr[12]
        rd['surface.k25rev'] = kfbnd * 59.0
        rd['surface.k27rev'] = rxntfr[13]
        rd['surface.k29rev'] = kfbnd * 0.1
        rd['surface.k31rev'] = rxntfr[14]
        rd['surface.k32rev'] = kfbnd * 1.0
        rd['surface.k33rev'] = rxntfr[15]
        rd['surface.k34rev'] = kfbnd * 0.07
        rd['surface.k35rev'] = rxntfr[16]

        for ii in ('1', '2', '4', '5', '10', '11', '13', '14', '16', '17', '22', '23', '25', '27', '29', '31', '32', '33', '34', '35'):
            rd['endosome.k' + ii + 'rev'] = rd['surface.k' + ii + 'rev'] * 5.0

        rd['endo'], rd['activeEndo'], rd['sortF'], rd['kRec'], rd['kDeg'] = tuple(rxntfr[17:22])
        rd['Rexpr_2Ra'], rd['Rexpr_2Rb'], rd['Rexpr_gc'], rd['Rexpr_15Ra'], rd['Rexpr_7R'], rd['Rexpr_9R'], rd['Rexpr_4Ra'], rd['Rexpr_21Ra'] = tuple(rxntfr[22:30])

    return rd


def getRateVec(rxntfr):
    """Retrieves and unpacks ordered dict + constructs rate vector for model fitting"""
    entries = rxntfr.size
    rxnlength = rxntfr.shape[0]

    if (entries / rxnlength) > 1:
        FullRateVec = np.zeros([rxntfr.shape[0], 60])
        for row in range(rxntfr.shape[0]):
            ratesParamsDict = getparamsdict(rxntfr[row, :])
            FullRateVec[row, :] = np.array(list(ratesParamsDict.values()), dtype=np.float)
    else:
        FullRateVec = np.zeros(60)
        ratesParamsDict = getparamsdict(rxntfr)
        FullRateVec = np.array(list(ratesParamsDict.values()), dtype=np.float)

    return FullRateVec
