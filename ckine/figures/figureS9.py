import os
import pandas as pd
import seaborn as sns
from os.path import join
from .figureCommon import subplotLabel, getSetup
from ..flow_meyer import process_sample_ILCs, form_gate
from FlowCytometryTools import FCMeasurement
from matplotlib import pyplot as plt
from matplotlib import cm

path_here = os.path.dirname(os.path.dirname(__file__))
plt.rcParams['svg.fonttype'] = 'none'


def makeFigure():
    """Get a list of the axis objects and create a figure"""

    ax, f = getSetup((7.5, 4.5), (2, 4))
    subplotLabel(ax)

    donor = "X"

    gateDF = pd.read_csv(join(path_here, "data/Meyer_Flow_Gates_ILCs.csv"))
    gateDF = gateDF.loc[gateDF.Donor == donor]
    sample = FCMeasurement(ID="X", datafile=(join(path_here, "data/Flow_Data_Meyer/ILCs/Donor_X/Specimen_001_Donor_X.fcs")))
    sample, _ = process_sample_ILCs(sample)
    sample = sample.subsample(0.3)

    PBMC_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "PBMC Gate")].Gate.values[0])
    Live_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "Live Gate")].Gate.values[0])
    T_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "T Gate")].Gate.values[0])
    ILC_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "ILC Gate")].Gate.values[0])
    Treg_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "Treg Gate")].Gate.values[0])
    ILC2_Gate = form_gate(gateDF.loc[(gateDF["Gate Label"] == "ILC2 Gate")].Gate.values[0])

    sample.plot(['FSC-H', 'SSC-H'], gates=PBMC_Gate, ax=ax[0], cmap=cm.jet, gate_colors=['red'])
    sample = sample.gate(PBMC_Gate)
    ax[0].set(title="Singlet Lymphocytes")

    sample.plot(['SSC-A', 'SSC-H'], gates=Live_Gate, ax=ax[1], cmap=cm.jet, gate_colors=['red'])
    sample = sample.gate(Live_Gate)
    ax[1].set(title="Live Gating")

    sample.plot(['FSC-A', 'Lin'], gates=T_Gate & ILC_Gate, ax=ax[2], cmap=cm.jet, gate_colors=['red', 'red'])
    ax[2].set(xlim=(60000, 150000), title="ILC Gating")
    T_sample = sample.gate(T_Gate)
    ILC_sample = sample.gate(ILC_Gate)

    T_sample.plot(['CD25', 'FoxP3'], gates=Treg_Gate, ax=ax[3], cmap=cm.jet, gate_colors=['red'])
    ax[3].set(title="Treg Gating")

    ILC_sample.data = ILC_sample.data.clip(lower=0, upper=10000)
    ILC_sample.plot(['CRTH2', 'CD127'], gates=ILC2_Gate, ax=ax[4], cmap=cm.jet, gate_colors=['red'])
    ax[4].set(ylim=(0, 10000), title="ILC Gating")
    
    SC_Data = pd.read_csv(join(path_here, "data/Meyer_Flow_ILCs_SC.csv"))
    sns.histplot(data=SC_Data[SC_Data.Cell == "Treg"], x="CD25", hue="Cell", ax=ax[5])
    ax[5].set(xlim=(0, 10000))
    sns.histplot(data=SC_Data[SC_Data.Cell == "ILC2"], x="CD25", hue="Cell", ax=ax[6])
    ax[6].set(xlim=(0, 10000))

    ILC_Data = pd.read_csv(join(path_here, "data/Meyer_Flow_ILCs.csv"))
    sns.barplot(data=ILC_Data, x="Donor", y="CD25", hue="Cell", ax=ax[7])
    ax[7].set(ylim=(0, 8000))
    ax[7].set_xticklabels(ax[7].get_xticklabels(), rotation=45, ha="right")
    plt.grid()

    return f
