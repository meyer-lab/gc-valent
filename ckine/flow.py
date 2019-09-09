"""
This file includes various methods for flow cytometry analysis.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from FlowCytometryTools import FCMeasurement
from FlowCytometryTools import QuadGate, ThresholdGate
from sklearn import preprocessing
from sklearn.decomposition import PCA


def importF(pathname):
    """
    Import FCS files. Variable input: name of path name to file. Output is a list of Data File Names in FCT Format
    Title/file names are returned in the array file --> later referenced in other functions as title/titles input argument
    """
    # Declare arrays and int
    file = []
    sample = []
    z = 0
    # Read in user input for file path and assign to array file
    pathlist = Path(r"" + str(pathname)).glob("**/*.fcs")
    for path in pathlist:
        path_in_str = str(path)
        file.append(path_in_str)
    file.sort()
    # Go through each file and assign the file contents to entry in the array sample
    for entry in file:
        sample.append(FCMeasurement(ID="Test Sample" + str(z), datafile=entry))
        z += 1
    # Returns the array sample which contains data of each file in folder (one file per entry in array)
    return sample, file


# *********************************** Gating Fxns *******************************************
# Treg and NonTreg

# add which channels relate to the proteins
def cd4():
    """Function for gating CD4+ cells (generates T cells)"""
    cd41 = ThresholdGate(6.514e+03, ('VL4-H'), region="above", name='cd41')
    cd42 = ThresholdGate(7.646e+03, ('VL4-H'), region="below", name='cd42')
    cd4_gate = cd41 & cd42
    return cd4_gate


def treg():
    """Function for creating and returning the T reg gate on CD4+ cells"""
    treg1 = QuadGate((4.814e+03, 3.229e+03), ('BL1-H', 'VL1-H'), region='top right', name='treg1')
    treg2 = QuadGate((6.258e+03, 5.814e+03), ('BL1-H', 'VL1-H'), region='bottom left', name='treg2')
    treg_gate = treg1 & treg2 & cd4()
    return treg_gate


def nonTreg():
    """Function for creating and returning the non T reg gate on CD4+ cells"""
    nontreg1 = QuadGate((5.115e+03, 3.470e+02), ('BL1-H', 'VL1-H'), region="top left", name='nontreg1')
    nontreg2 = QuadGate((2.586e+03, 5.245e+03), ('BL1-H', 'VL1-H'), region="bottom right", name='nontreg2')
    nonTreg_gate = nontreg1 & nontreg2 & cd4()
    return nonTreg_gate


def nk():
    """Function for creating and returning the NK gate"""
    # NK cells: Take quad gates for NK cells and combine them to create single, overall NK gate
    nk1 = QuadGate((6.468e03, 4.861e03), ("BL1-H", "VL4-H"), region="top left", name="nk1")
    nk2 = QuadGate((5.550e03, 5.813e03), ("BL1-H", "VL4-H"), region="bottom right", name="nk2")
    nk_gate = nk1 & nk2
    return nk_gate


def bnk():
    """Function for creating and returning the BNK gate"""
    # Bright NK cells: Take quad gates for bright NK cells and combine them to create single, overall bright NK gate
    bnk1 = QuadGate((7.342e03, 4.899e03), ("BL1-H", "VL4-H"), region="top left", name="bnk1")
    bnk2 = QuadGate((6.533e03, 5.751e03), ("BL1-H", "VL4-H"), region="bottom right", name="bnk2")
    bnk_gate = bnk1 & bnk2
    return bnk_gate


def cd():
    """Function for creating and returning the CD gate"""
    # CD cells: Take quad gates for CD cells and combine them to create single, overall CD gate
    cd1 = QuadGate((9.016e03, 5.976e03), ("RL1-H", "VL4-H"), region="top left", name="cd1")
    cd2 = QuadGate((6.825e03, 7.541e03), ("RL1-H", "VL4-H"), region="bottom right", name="cd2")
    cd_gate = cd1 & cd2
    return cd_gate


def cellCount(sample_i, gate):
    """
    Function for returning the count of cells in a single .fcs. file of a single cell file. Arguments: single sample/.fcs file and the gate of the
    desired cell output.
    """
    # Import single file and save data to a variable --> transform to logarithmic scale
    smpl = sample_i.transform("hlog", channels=["BL1-H", "VL1-H", "VL4-H", "RL1-H"])
    # Apply T reg gate to overall data --> i.e. step that detrmines which cells are T reg
    cells = smpl.gate(gate)
    # Number of events (AKA number of cells)
    cell_count = cells.get_data().shape[0]
    # print('Number of Treg cells:' + str(treg_count))
    return cell_count


def rawData(sample_i, gate):
    """
    Function that returns the raw data of certain cell population in a given file. Arguments: sample_i is a single entry/.fcs file and the gate
    of the desired cell population.
    """
    smpl = sample_i.transform("hlog", channels=["BL1-H", "VL1-H", "VL4-H", "RL1-H"])
    # Apply T reg gate to overall data --> i.e. step that detrmines which cells are T reg
    cells = smpl.gate(gate)
    # Get raw data of t reg cells in file
    cell_data = cells.get_data()
    return cell_data


def tcells(sample_i, treg_gate, nonTreg_gate, title):
    """
    Function that is used to plot the Treg and NonTreg gates in CD4+ cells. Treg (yellow) and Non Treg (green). sample_i is an indivual flow cytommetry file/data.
    """
    # Data to use is on CD4+ cells
    # Apply new T reg and Non treg gate
    # Assign data of current file for analysis to variable smpl and transform to log scale
    smpl = sample_i.transform('hlog', channels=["VL4-H", "BL1-H", "VL1-H"])
    # Create data set to only include CD4 cells
    cd4_gate = cd4()
    cd4_cells = smpl.gate(cd4_gate)
    # CD25 v. Foxp33: VL1 v. BL1
    # Treg
    # Apply T reg gate to overall data --> step that determines which cells are Treg
    treg_cells = smpl.gate(treg_gate)
    # Non Tregs
    # Apply non T reg gate to overall data --> step that detrmines which cells are non T reg
    nonTreg_cells = smpl.gate(nonTreg_gate)

    # Declare figure and axis
    _, ax = plt.subplots()
    # Plot the treg gate
    treg_cells.plot(["BL1-H", "VL1-H"], color="teal")
    # Plot the non Treg gate
    nonTreg_cells.plot(["BL1-H", "VL1-H"], color="cyan")
    # Plot all of the cells in the file
    ax.set_title("T Reg + Non T Reg - Gating - " + str(title), fontsize=12)
    cd4_cells.plot(["BL1-H", "VL1-H"])
    plt.xlabel("Foxp3", fontsize=12)
    plt.ylabel("CD25", fontsize=12)
    # Set values for legend
    bar_T = ax.bar(np.arange(0, 10), np.arange(1, 11), color="teal")
    bar_NT = ax.bar(np.arange(0, 10), np.arange(30, 40), bottom=np.arange(1, 11), color="cyan")
    ax.legend([bar_T, bar_NT], ("T Reg", "Non T Reg"), loc="upper left")


def nk_bnk_plot(sample_i, nk_gate, bnk_gate, title):
    """
    Function that plots the graph of NK and Bright NK cells (both are determined by same x, y-axis). Arguemnt 1: current sample (a single file).
    Argument 2: the gate for NK. Argument 3: the gate for bright NK.
    """
    smpl = sample_i.transform("hlog", channels=["BL1-H", "VL4-H", "RL1-H"])

    # CD3 v. CD56: VL4 v. BL1
    # NK
    # Apply NK gate to overall data --> step that determines which cells are NK
    nk_cells = smpl.gate(nk_gate)
    # CD56 Bright NK
    # Apply Bright NK gate to overall data --> step that determines which cells are Bright NK
    bnk_cells = smpl.gate(bnk_gate)

    _, ax1 = plt.subplots()
    ax1.set_title("CD56 BrightNK + NK - Gating - " + str(title), fontsize=12)
    nk_cells.plot(["BL1-H", "VL4-H"], color="y", label="NK")
    bnk_cells.plot(["BL1-H", "VL4-H"], color="g", label="Bright NK")
    smpl.plot(["BL1-H", "VL4-H"])

    bar_NK = ax1.bar(np.arange(0, 10), np.arange(1, 11), color="y")
    bar_BNK = ax1.bar(np.arange(0, 10), np.arange(30, 40), bottom=np.arange(1, 11), color="g")
    ax1.legend([bar_NK, bar_BNK], ("NK", "Bright NK"), loc="upper left")


def cd_plot(sample_i, cd_gate, title):
    """
    Function that plots the graph of CD cells. Argument 1: current sample (a single file). Argument 2: the gate for CD cells. Argument 3: the value
    of the current i in a for loop --> use
    when plotting multiple files.
    """
    smpl = sample_i.transform("hlog", channels=["BL1-H", "VL4-H", "RL1-H"])
    # CD3 v. CD8: VL4 v. RL1
    # CD3+CD8+
    # Apply CD cell gate to overall data --> step that determines which cells are CD
    cd_cells = smpl.gate(cd_gate)

    _, ax2 = plt.subplots()
    ax2.set_title("CD3+CD8+ - Gating - " + str(title), fontsize=20)
    cd_cells.plot(["RL1-H", "VL4-H"], color="b")
    smpl.plot(["RL1-H", "VL4-H"])

    bar_CD = ax2.bar(np.arange(0, 10), np.arange(1, 11), color="b")
    ax2.legend([bar_CD], ("CD3+8+"), loc="upper left")


def count_data(sampleType, gate):
    """
    Used to count the number of cells and store the data of all of these cells in a folder with multiple files --> automates the process sampleType
    is NK or T cell data, gate is the desired cell population.
    Sample type: is the overall importF assignment for T or NK (all the T cell files, all NK cell files)
    """
    # declare the arrays to store the data
    count_array = []
    data_array = []
    # create the for loop to file through the data and save to the arrays
    # using the functions created above for a singular file
    for _, sample in enumerate(sampleType):
        count_array.append(cellCount(sample, gate))
        data_array.append(rawData(sample, gate))
    # returns the array for count of cells and the array where each entry is the data for the specific cell population in that .fcs file
    return count_array, data_array


def plotAll(sampleType, check, gate1, gate2, titles):
    """
    Ask the user to input 't' for t cell, 'n' for nk cell, and 'c' for cd cell checks are used to determine if user input a T-cell, NK-cell, or
    CD-cell gate automates the process for plotting multiple files.
    """
    if check == "t":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            tcells(sample, gate1, gate2, title)
    elif check == "n":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            nk_bnk_plot(sample, gate1, gate2, title)
    elif check == "c":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            cd_plot(sample, gate1, title)


# ********************************** PCA Functions****************************************************
def sampleT(smpl):
    """Output is the T cells data (the protein channels related to T cells)"""
    # Features are the protein channels of interest when analyzing T cells
    features = ["BL1-H", "VL1-H", "VL4-H", "BL3-H"]
    # Transform to put on log scale
    tform = smpl.transform("hlog", channels=["BL1-H", "VL1-H", "VL4-H", "BL3-H", "RL1-H"])
    # Save the data of each column of the protein channels
    data = tform.data[["BL1-H", "VL1-H", "VL4-H", "BL3-H"]][0:]
    # Save pSTAT5 data
    pstat = tform.data[["RL1-H"]][0:]
    return data, pstat, features


def sampleNK(smpl):
    """Output is the NK cells data (the protein channels related to NK cells)"""
    # For NK, the data consists of different channels so the data var. output will be different
    # Output is data specific to NK cells
    # Features for the NK file of proteins (CD3, CD8, CD56)
    features = ["VL4-H", "RL1-H", "BL1-H"]
    # Transform all proteins (including pSTAT5)
    tform = smpl.transform("hlog", channels=["VL4-H", "RL1-H", "BL1-H", "BL2-H"])
    # Assign data of three protein channels AND pSTAT5
    data = tform.data[["VL4-H", "RL1-H", "BL1-H"]][0:]
    pstat = tform.data[["BL2-H"]][0:]
    return data, pstat, features


def fitPCA(data, features):
    """
    Fits the PCA model to data, and returns the fit PCA model for future transformations
    (allows for consistent loadings plots between wells)
    """
    # Apply PCA to the data set
    # setting values of data of selected features to data frame
    xi = data.loc[:, features].values
    # STANDARDIZE DATA --> very important to do before applying machine learning algorithm
    scaler = preprocessing.StandardScaler()
    xs = scaler.fit_transform(xi)
    xs = np.nan_to_num(xs)
    # setting how many components wanted --> PC1 and PC2
    pca = PCA(n_components=2)
    # apply PCA to standardized data set
    PCAobj = pca.fit(xs)
    # creates the loading array (equation is defintion of loading)
    loading = pca.components_.T
    return PCAobj, loading


def appPCA(data, features, PCAobj):
    """Applies the PCA algorithm to the data set"""
    # setting values of data of selected features to data frame
    xi = data.loc[:, features].values
    # STANDARDIZE DATA --> very important to do before applying machine learning algorithm
    scaler = preprocessing.StandardScaler()
    xs = scaler.fit_transform(xi)
    xs = np.nan_to_num(xs)
    # transform to the prefit pca object
    xf = PCAobj.transform(xs)
    return xf


def pcaPlt(xf, pstat, features, title, tplate=True):
    """
    Used to plot the score graph.
    Scattered point color gradients are based on range/abundance of pSTAT5 data. Light --> Dark = Less --> More Active
    """
    # PCA
    if len(features) == 4:
        name = "T Cells"
    elif len(features) == 3:
        name = "NK Cells"
    # Setting x and y values from xf
    x = xf[:, 0]
    y = xf[:, 1]
    # Saving numerical values for pSTAT5 data
    pstat_data = pstat.values
    # Creating a figure for both scatter and mesh plots for PCA
    # This is the scatter plot of the cell clusters colored by pSTAT5 data
    # lighter --> darker = less --> more pSTAT5 present
    # Creating correct dimensions
    pstat = np.squeeze(pstat_data)
    # Determining overall pSTAT5 for a certain fcs data set
    pstat_mean = np.mean(pstat)
    print("mean pStat activity: ", pstat_mean)
    # Creating a data from of x, y, and pSTAT5 in order to graph using seaborn
    combined = np.stack((x, y, pstat)).T
    df = pd.DataFrame(combined, columns=["PC1", "PC2", "pSTAT5"])
    # Creating plot using seaborn. Cool note: virdis is visible for individuals who are colorblind.
    _, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(name + " - PCA - " + str(title), fontsize=20)
    plt.xlim(-4, 6)
    plt.ylim(-4, 4)
    if tplate:
        sns.scatterplot(x="PC1", y="PC2", hue="pSTAT5", palette="viridis", data=df, s=10, ax=ax, legend=False, hue_norm=(3000, 7000))
        points = plt.scatter(df["PC1"], df["PC2"], c=df["pSTAT5"], s=0, cmap="viridis", vmin=3000, vmax=7000) #set style options
    else:
        sns.scatterplot(x="PC1", y="PC2", hue="pSTAT5", palette="viridis", data=df, s=10, ax=ax, legend=False, hue_norm=(0, 5000))
        points = plt.scatter(df["PC1"], df["PC2"], c=df["pSTAT5"], s=0, cmap="viridis", vmin=0, vmax=5000) #set style options
    ax.set_xlabel("PC1", fontsize=15)
    ax.set_ylabel("PC2", fontsize=15)
    #add a color bar
    plt.colorbar(points)


def loadingPlot(loading, features, i, title):
    """Plot the loading data"""
    # Loading
    # Create graph for loading values
    x_load = loading[:, 0]
    y_load = loading[:, 1]

    # Create figure for the loading plot
    fig1 = plt.figure(figsize=(8, 8))
    ax = fig1.add_subplot(1, 1, 1)
    ax.set_xlabel("PC1", fontsize=15)
    ax.set_ylabel("PC2", fontsize=15)
    plt.scatter(x_load, y_load)

    for z, feature in enumerate(features):
        # Please note: not the best logic, but there are three features in NK and four features in T cells
        if len(features) == 4:
            name = "T Cells"
            if feature == "BL1-H":
                feature = "Foxp3"
            elif feature == "VL1-H":
                feature = "CD25"
            elif feature == "VL4-H":
                feature = "CD4"
            elif feature == "BL3-H":
                feature = "CD45RA"
        if len(features) == 3:
            name = "NK Cells"
            if feature == "VL4-H":
                feature = "CD3"
            if feature == "RL1-H":
                feature = "CD8"
            if feature == "BL1-H":
                feature = "CD56"
        plt.annotate(str(feature), xy=(x_load[z], y_load[z]))
        plt.savefig('loading' + str(i) + '.png')
    ax.set_title(name + " - Loading - " + str(title), fontsize=20)


def pcaAll(sampleType, check, titles):
    """
    Use to plot the score and loading graphs for PCA. Assign protein and pstat5 arrays AND score and loading arrays
    This is all the data for each file.
    Want to use for both T and NK cells? Use it twice!
    sampleType is importF for T or NK
    check == "t" for T cells OR check == "n" for NK cells
    """
    # declare the arrays to store the data
    data_array = []
    pstat_array = []
    xf_array = []
    # create the for loop to file through the data and save to the arrays
    # using the functions created above for a singular file
    if check == "t":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            data, pstat, features = sampleT(sample)
            data_array.append(data)
            pstat_array.append(pstat)
            if i == 0:
                PCAobj, loading = fitPCA(data, features)
            xf = appPCA(data, features, PCAobj)
            xf_array.append(xf)
            pcaPlt(xf, pstat, features, title, tplate=True)
            loadingPlot(loading, features, i, title)
            plt.show()
    elif check == "n":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            data, pstat, features = sampleNK(sample)
            data_array.append(data)
            pstat_array.append(pstat)
            if i == 0:
                PCAobj, loading = fitPCA(data, features)
            xf = appPCA(data, features, PCAobj)
            pcaPlt(xf, pstat, features, title, tplate=False)
            loadingPlot(loading, features, i, title)
            plt.show()
    return data_array, pstat_array, xf_array

# ************************PCA by color (gating+PCA)******************************


def sampleTcolor(smpl):
    """Output is the T cells data (the protein channels related to T cells)"""
    # Features are the protein channels of interest when analyzing T cells
    features = ["BL1-H", "VL1-H", "VL4-H", "BL3-H"]
    # Transform to put on log scale
    cd45dat = smpl[["BL3-H"]]
    cd45dat = cd45dat.iloc[:, 0]
    cd45dat = np.log10(cd45dat)
    tform = smpl.transform("hlog", channels=["BL1-H", "VL1-H", "VL4-H", "BL3-H", "RL1-H"])
    # Save the data of each column of the protein channels
    data = tform.data[["BL1-H", "VL1-H", "VL4-H", "BL3-H"]][0:]
    # Save pSTAT5 data
    pstat = tform.data[["RL1-H"]][0:]
    colmat = []*(len(data)+1)
    for i in range(len(data)):
        if data.iat[i, 0] > 4.814e+03 and data.iat[i, 0] < 6.258e+03 and data.iat[i, 1] > 3.229e+03 and data.iat[i, 1] < 5.814e+03:
            if cd45dat[i] > 5:
                colmat.append('r') #Treg naive
            else:
                colmat.append('darkorange') #Treg mem
        elif data.iat[i, 0] > 2.586e+03 and data.iat[i, 0] < 5.115e+03 and data.iat[i, 1] > 3.470e+02 and data.iat[i, 1] < 5.245e+03:
            if cd45dat[i] > 5:
                colmat.append('g') #Thelp naive
            else:
                colmat.append('darkorchid') #Thelp mem
        else:
            colmat.append('c')
    return data, pstat, features, colmat


def sampleNKcolor(smpl):
    """Output is the NK cells data (the protein channels related to NK cells)"""
    # For NK, the data consists of different channels so the data var. output will be different
    # Output is data specific to NK cells
    # Features for the NK file of proteins (CD3, CD8, CD56)
    features = ["VL4-H", "RL1-H", "BL1-H"]
    # Transform all proteins (including pSTAT5)
    tform = smpl.transform("hlog", channels=["VL4-H", "RL1-H", "BL1-H", "BL2-H"])
    # Assign data of three protein channels AND pSTAT5
    data = tform.data[["VL4-H", "RL1-H", "BL1-H"]][0:]
    pstat = tform.data[["BL2-H"]][0:]
    # Create a section for assigning colors to each data point of each cell population --> in this case NK cells
    colmat = [] * (len(data) + 1)
    for i in range(len(data)):
        if data.iat[i, 0] > 5.550e03 and data.iat[i, 0] < 6.468e03 and data.iat[i, 2] > 4.861e03 and data.iat[i, 2] < 5.813e03:
            colmat.append('r')  # nk
        elif data.iat[i, 0] > 6.533e03 and data.iat[i, 0] < 7.34e03 and data.iat[i, 2] > 4.899e03 and data.iat[i, 2] < 5.751e03:
            colmat.append('g')  # bnk
        else:
            colmat.append('c')
    return data, pstat, features, colmat


def pcaPltColor(xf, pstat, features, title, colormat, tregc=True):
    """
    Used to plot the score graph.
    Scattered point color gradients are based on range/abundance of pSTAT5 data. Light --> Dark = Less --> More Active
    """
    # PCA
    if len(features) == 4:
        name = "T Cells"
    elif len(features) == 3:
        name = "NK Cells"
    # Setting x and y values from xf
    x = xf[:, 0]
    y = xf[:, 1]
    # Working with pSTAT5 data --> setting min and max values
    pstat_data = pstat.values
    pstat_mean = np.mean(pstat_data)
    print("mean pStat activity: ", pstat_mean)
    plt.show()
    # Creating a figure for both scatter and mesh plots for PCA
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Principal Component 1", fontsize=15)
    ax.set_ylabel("Principal Component 2", fontsize=15)
    ax.set_title(name + " - PCA - " + str(title), fontsize=20)
    ax.set(xlim=(-5, 5), ylim=(-5, 5))
    # This is the scatter plot of the cell clusters colored by pSTAT5 data
    # lighter --> darker = less --> more pSTAT5 present
    colormat = np.array(colormat)
    if tregc:
        plt.scatter(x[colormat == "c"], y[colormat == "c"], s=1, c="c", label="Other", alpha=0.5)
        plt.scatter(x[colormat == "g"], y[colormat == "g"], s=1, c="g", label="T Helper Naive", alpha=0.5)
        plt.scatter(x[colormat == "darkorchid"], y[colormat == "darkorchid"], s=1, c="darkorchid", label="T Helper Memory", alpha=0.5)
        plt.scatter(x[colormat == "darkorange"], y[colormat == "darkorange"], s=1, c="darkorange", label="T Reg Memory", alpha=0.5)
        plt.scatter(x[colormat == "r"], y[colormat == "r"], s=1, c="r", label="T Reg Naive", alpha=0.5)
        plt.legend()
    else:
        plt.scatter(x[colormat == "c"], y[colormat == "c"], s=1, c="c", label="Other", alpha=0.5)
        plt.scatter(x[colormat == "g"], y[colormat == "g"], s=1, c="g", label="BNK", alpha=0.5)
        plt.scatter(x[colormat == "r"], y[colormat == "r"], s=1, c="r", label="NK", alpha=0.5)
        plt.legend()


def pcaAllCellType(sampleType, check, titles):
    """
    Use to plot the score and loading graphs for PCA. Assign protein and pstat5 arrays AND score and loading arrays
    This is all the data for each file.
    Want to use for both T and NK cells? Use it twice!
    sampleType is importF for T or NK
    check == "t" for T cells OR check == "n" for NK cells
    """
    # declare the arrays to store the data
    data_array = []
    pstat_array = []
    xf_array = []
    loading_array = []

    # create the for loop to file through the data and save to the arrays
    # using the functions created above for a singular file
    if check == "t":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            data, pstat, features, colormat = sampleTcolor(sample)
            data_array.append(data)
            pstat_array.append(pstat)
            if i == 0:
                PCAobj, loading = fitPCA(data, features)
            xf = appPCA(data, features, PCAobj)
            xf_array.append(xf)
            loading_array.append(loading)
            pcaPltColor(xf, pstat, features, title, colormat)  # changed
            loadingPlot(loading, features, i, title)
    elif check == "n":
        for i, sample in enumerate(sampleType):
            title = titles[i].split("/")
            title = title[len(title) - 1]
            data, pstat, features, colormat = sampleNKcolor(sample)
            data_array.append(data)
            pstat_array.append(pstat)
            if i == 0:
                PCAobj, loading = fitPCA(data, features)
            xf = appPCA(data, features, PCAobj)
            pcaPltColor(xf, pstat, features, title, colormat)
            loadingPlot(loading, features, i, title)
    plt.show()
    return data_array, pstat_array, xf_array, loading_array

#************************Dose Response by PCA******************************

def PCADoseResponse(sampleType, PC1Bnds, PC2Bnds, Timepoint):
    """ Given data from a time Point and two PC bounds, the dose response curve will be calculated and graphed (needs folder with FCS from one time point)"""
    dosemat = np.array([84, 28, 9.333333, 3.111, 1.037037, 0.345679, 0.115226, 0.038409, 0.012803, 0.004268, 0.001423, 0.000474])
    Pstatvals = []

    for i, sample in enumerate(sampleType):
        data, pstat, features = sampleT(sample) #retrieve data
        if i == 0:
            PCAobj, _ = fitPCA(data, features) #only fit to first set
        xf = appPCA(data, features, PCAobj) #get PC1/2 vals
        PC1, PC2, pstat = np.transpose(xf[:, 0]), np.transpose(xf[:, 1]), pstat.to_numpy()
        PC1, PC2 = np.reshape(PC1, (PC1.size, 1)), np.reshape(PC2, (PC2.size, 1))
        PCAstat = np.concatenate((PC1, PC2, pstat), axis=1)
        PCApd = pd.DataFrame({'PC1': PCAstat[:, 0], 'PC2': PCAstat[:, 1], 'Pstat': PCAstat[:, 2]}) #arrange into pandas datafrome
        PCApd = PCApd[PCApd['PC1'] >= PC1Bnds[0]] #remove data that that is not within given PC bounds
        PCApd = PCApd[PCApd['PC1'] <= PC1Bnds[1]]
        PCApd = PCApd[PCApd['PC2'] >= PC2Bnds[0]]
        PCApd = PCApd[PCApd['PC2'] <= PC2Bnds[1]]
        Pstatvals.append(PCApd.loc[:, "Pstat"].mean()) #take average Pstat activity of data fitting criteria

    _, ax = plt.subplots(figsize=(8, 8))
    plt.scatter(dosemat, Pstatvals)
    ax.set_title(Timepoint + " PCA Gated Dose Response Curve", fontsize=20)
    ax.set_xscale('log')
    ax.set_xlabel("Cytokine Dosage (log10[nM])", fontsize=15)
    ax.set_ylabel("Average Pstat Activity", fontsize=15)
    ax.set(xlim=(0.0001, 100))
    plt.show()
