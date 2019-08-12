"""
This file includes various methods for flow cytometry analysis.
"""
#import FlowCytometryTools

#*********************************** Gating Fxns *******************************************
# Import FCS files. Variable input: name of path name to file. 
def importF(pathname):
    # Declare arrays and int
    file = []
    sample = []
    z = 0
    # Read in user input for file path and assign to array file 
    pathlist = Path(r'' + str(pathname)).glob('**/*.fcs')
    for path in pathlist:
        path_in_str = str(path)
        file.append(path_in_str)
    # Go through each file and assign the file contents to entry in the array sample
    for entry in file:
        sample.append(FCMeasurement(ID = 'Test Sample' + str(z), datafile = entry))
        z+=1
    importF.sample = sample
    # Returns the array sample which contains data of each file in folder (one file per entry in array)
    return sample
# Treg and NonTreg

# Function for creating and returning the T reg gate
def treg():
    # T reg: Take quad gates for T reg cells and combine them to create single, overall T reg gate
    treg1 = QuadGate((7.063e+03, 3.937e+03), ('BL1-H', 'VL1-H'), region='top left', name='treg1')
    treg2 = QuadGate((5.412e+03, 6.382e+03), ('BL1-H', 'VL1-H'), region='bottom right', name='treg2')
    treg_gate = treg1 & treg2
    return treg_gate

# Function for creating and returning the Non T reg gate
def nonTreg():
    # non T reg: Take quad gates for non T reg cells and combine to create overall non T reg gate
    nonTreg1 = QuadGate((5.233e+03, 1.731e+03), ('BL1-H', 'VL1-H'), region='top left', name='nonTreg1')
    nonTreg2 = QuadGate((2.668e+03, 5.692e+03), ('BL1-H', 'VL1-H'), region='bottom right', name='nonTreg2')
    nonTreg_gate = nonTreg1 & nonTreg2
    return nonTreg_gate

def nk():
    # NK cells: Take quad gates for NK cells and combine them to create single, overall NK gate
    nk1 = QuadGate((6.468e+03, 4.861e+03), ('BL1-H', 'VL4-H'), region='top left', name='nk1')
    nk2 = QuadGate((5.550e+03, 5.813e+03), ('BL1-H', 'VL4-H'), region='bottom right', name='nk2')
    nk_gate = nk1 & nk2
    return nk_gate

def bnk():
    # Bright NK cells: Take quad gates for bright NK cells and combine them to create single, overall bright NK gate
    bnk1 = QuadGate((7.342e+03, 4.899e+03), ('BL1-H', 'VL4-H'), region='top left', name='bnk1')
    bnk2 = QuadGate((6.533e+03, 5.751e+03), ('BL1-H', 'VL4-H'), region='bottom right', name='bnk2')
    bnk_gate = bnk1 & bnk2
    return bnk_gate

def cd():
    # CD cells: Take quad gates for CD cells and combine them to create single, overall CD gate
    cd1 = QuadGate((9.016e+03, 5.976e+03), ('RL1-H', 'VL4-H'), region='top left', name='cd1')
    cd2 = QuadGate((6.825e+03, 7.541e+03), ('RL1-H', 'VL4-H'), region='bottom right', name='cd2')
    cd_gate = cd1 & cd2
    return cd_gate

def cellCount(sample_i, gate):
    # Import single file and save data to a variable --> transform to logarithmic scale
    smpl = sample_i.transform('hlog', channels=['BL1-H','VL1-H','VL4-H','RL1-H'])
    # Apply T reg gate to overall data --> step that detrmines which cells are T reg
    cells = smpl.gate(gate)
    # Number of events (AKA number of cells)
    cell_count = cells.get_data().shape[0]
    #print('Number of Treg cells:' + str(treg_count))
    return cell_count

# Function that returns the raw data of certain cell population in a given file
def rawData(sample_i, gate):
    smpl = sample_i.transform('hlog', channels=['BL1-H','VL1-H','VL4-H','RL1-H'])
    # Apply T reg gate to overall data --> step that detrmines which cells are T reg
    cells = smpl.gate(gate)  
    # Get raw data of t reg cells in file
    cell_data = cells.get_data()
    return cell_data

# Function that is used to plot the Treg and NonTreg gates
# Treg (yellow) and Non Treg (green)
# sample_i is an indivual flow cytommetry file/data
def tcells(sample_i, treg_gate, nonTreg_gate):
    # Assign data of current file for analysis to variable smpl and transform to log scale
    smpl = sample_i.transform('hlog', channels=['BL1-H','VL1-H'])
    # CD25 v. Foxp33: VL1 v. BL1
    # Treg
    # Apply T reg gate to overall data --> step that determines which cells are Treg
    treg = smpl.gate(treg_gate)
    # Non Tregs
    # Apply non T reg gate to overall data --> step that detrmines which cells are non T reg
    nonTreg = smpl.gate(nonTreg_gate)
    
    # Declare figure and axis
    fig, ax = plt.subplots()
    # Plot the treg gate
    treg.plot(['BL1-H','VL1-H'], color = 'y')
    # Plot the non Treg gate
    nonTreg.plot(['BL1-H','VL1-H'], color = 'g')
    # Plot all of the cells in the file
    ax.set_title('T Reg + Non T Reg - Gating - File ' + str(i), fontsize = 20)
    smpl.plot(['BL1-H','VL1-H'])
      
    bar_T = ax.bar(np.arange(0,10), np.arange(1,11), color="y")
    bar_NT = ax.bar(np.arange(0,10), np.arange(30,40), bottom=np.arange(1,11), color="g")
    ax.legend([bar_T, bar_NT], ("T Reg", "Non T Reg"), loc='upper left')
    
# Function that plots the graph of NK and Bright NK cells (both are determined by same x, y-axis)
# Arguemnt 1: current sample (a single file)
# Argument 2: the gate for NK
# Argument 3: the gate for bright NK 

#combine to use NKcells, use cd_gate as argument
def nk_bnk_plot(sample_i, nk_gate, bnk_gate):
    smpl = sample_i.transform('hlog', channels=['BL1-H','VL4-H','RL1-H'])

    # CD3 v. CD56: VL4 v. BL1
    # NK
    # Apply NK gate to overall data --> step that determines which cells are NK
    nk = smpl.gate(nk_gate)
    # CD56 Bright NK
    # Apply Bright NK gate to overall data --> step that determines which cells are Bright NK
    bnk = smpl.gate(bnk_gate)
    
    fig1, ax1 = plt.subplots()
    ax1.set_title('CD56 BrightNK + NK - Gating - File ' + str(i), fontsize = 20)
    nk.plot(['BL1-H','VL4-H'], color = 'y', label = 'NK')
    bnk.plot(['BL1-H','VL4-H'], color = 'g', label = 'Bright NK')
    smpl.plot(['BL1-H','VL4-H'])
    
    bar_NK = ax1.bar(np.arange(0,10), np.arange(1,11), color="y")
    bar_BNK = ax1.bar(np.arange(0,10), np.arange(30,40), bottom=np.arange(1,11), color="g")
    ax1.legend([bar_NK, bar_BNK], ("NK", "Bright NK"), loc='upper left')
    
# Function that plots the graph of CD cells
# Argument 1: current sample (a single file)
# Argument 2: the gate for CD cells
def cd_plot(sample_i, cd_gate):
    smpl = sample_i.transform('hlog', channels=['BL1-H','VL4-H','RL1-H'])
    #CD3 v. CD8: VL4 v. RL1
    #CD3+CD8+
    # Apply CD cell gate to overall data --> step that determines which cells are CD
    cd = smpl.gate(cd_gate)

    fig2, ax2 = plt.subplots()
    ax2.set_title('CD3+CD8+ - Gating - File ' + str(i), fontsize = 20)
    cd.plot(['RL1-H','VL4-H'],color='b')
    smpl.plot(['RL1-H','VL4-H'])
    
    bar_CD = ax2.bar(np.arange(0,10), np.arange(1,11), color='b')
    ax2.legend([bar_CD], ('CD3+8+'), loc='upper left')
    
#********************************** PCA Functions****************************************************
def sampleT(smpl): 
    
    features = ['BL1-H', 'VL1-H', 'VL4-H', 'BL3-H']
    # Output is data specific to T cells
    # Transform to put on log scale
    tform = smpl.transform('hlog', channels=['BL1-H', 'VL1-H', 'VL4-H', 'BL3-H', 'RL1-H'])
    # Save the data of each column of the protein channels
    data = tform.data[['BL1-H', 'VL1-H', 'VL4-H', 'BL3-H']][0:]
    # Save pSTAT5 data
    pstat = tform.data[['RL1-H']][0:]
    
    return data, pstat, features

def sampleNK(smpl):
    
    # For NK, the data consists of different channels so the data var. output will be different
    # Output is data specific to NK cells
    # Features for the NK file of proteins (CD3, CD8, CD56)
    features = ['VL4-H','RL1-H','BL1-H']
    # Transform all proteins (including pSTAT5)
    tform = smpl.transform('hlog', channels=['VL4-H','RL1-H','BL1-H','BL2-H'])
    # Assign data of three protein channels AND pSTAT5
    data = tform.data[['VL4-H','RL1-H','BL1-H']][0:]
    pstat = tform.data[['BL2-H']][0:]
    
    return data, pstat, features

def appPCA(data, features):
    
    # Apply PCA to the data set
    # setting values of data of selected features to data frame
    xi = data.loc[:,features].values
    # STANDARDIZE DATA --> very important to do before applying machine learning algorithm
    xs = sklearn.preprocessing.scale(xi)
    xs = np.nan_to_num(xs)
    # setting how many components wanted --> PC1 and PC2
    pca = PCA(n_components=2)
    # apply PCA to standardized data set
    # NOTE: score == xf
    xf = pca.fit(xs).transform(xs)
    
    loading = pca.components_.T * np.sqrt(pca.explained_variance_)

    # These are the values for each cell on our new PC1 and PC2 axis graph
    principalDf = pandas.DataFrame(data = xf, columns = ['principal component 1', 'principal component 2'])
    
    return xf, loading

def pcaPltCat(xf, pstat, loading, features):
    
    ### WORK on how to get the two scatter plots added onto the SAME figure (look up proper syntax for declaring figures)
    
    # Plot type #2
    # ************** K-means --> replace with other type of PCA
    
    # Create graph for loading values
    x_load = loading[:,0]
    y_load = loading[:,1]   
    
    # Create figure for the loading plot
    fig1 = plt.figure(figsize = (8,8))
    ax = fig1.add_subplot(1,1,1)
    ax.set_xlabel('PC1', fontsize = 15)
    ax.set_ylabel('PC2', fontsize = 15)
    #ax.set(xlim=(-1, 1), ylim=(-1,1))
    
    plt.scatter(x_load, y_load)
    
    for z, feature in enumerate(features):
        # Please note: not the best logic, but there are three features in NK and four features in T cells
        if len(features) == 4:
            name = 'T Cells'
            if feature == 'BL1-H':
                feature = 'Foxp3'
            elif feature == 'VL1-H':
                feature = 'CD25'
            elif feature == 'VL4-H':
                feature = 'CD4'
            elif feature == 'BL3-H':
                feature = 'CD45RA'       
        if len(features) == 3:
            name = 'NK Cells'
            if feature == 'VL4-H':
                feature = 'CD3'
            if feature == 'RL1-H':
                feature = 'CD8'
            if feature == 'BL1-H':
                feature = 'CD56'
        plt.annotate(str(feature), xy = (x_load[z], y_load[z]))
    ax.set_title(name + ' - Loading - File ' + str(i), fontsize = 20)
    
        
    #Categorize the different proteins -- four categories
    kmeans = sklearn.cluster.KMeans(n_clusters = 4)
    x_kmeans = kmeans.fit(xf)
    y_kmeans = kmeans.predict(xf)
    
    # Step size of the mesh. Decrease to increase the quality of the VQ.
    h = .02     
    # point in the mesh [x_min, x_max]x[y_min, y_max]
    x = xf[:,0]
    y = xf[:,1]
    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    # Create the code for the mesh grid background using the values of the boundaries
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    # Obtain labels for each point in mesh. Use last trained model.
    Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])
    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    
    # Working with pSTAT5 data --> setting min and max values
    p_min = pstat.values.min()
    p_max = pstat.values.max()
    pstat_data = pstat.values
    
    # Creating a figure for both scatter and mesh plots for PCA
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title(name + ' - PCA - File ' + str(i), fontsize = 20)
    ax.set(xlim=(-5, 5), ylim=(-5, 5))
    
    # This is the scatter plot of the cell clusters colored by pSTAT5 data
    # lighter --> darker = less --> more pSTAT5 present
    plt.scatter(x,y, s = 0.1, c=np.squeeze(pstat_data), cmap = "Greens")
    
    # Plot the centroids as a white X
    centroids = kmeans.cluster_centers_
    # Plot the center of the kmeans cluster ** will be replaced by PCA cluster/mesh
    plt.scatter(centroids[:, 0], centroids[:, 1], marker='x', s=169, linewidths=3, color='w', zorder=10) 
    # Displays the mesh category plot
    plt.imshow(Z, interpolation='nearest', extent=(xx.min(), xx.max(), yy.min(), yy.max()),cmap=plt.cm.Paired,aspect='auto', origin='lower')
    
# ********************** Gating Sample Program ***********************
# Import all necessary packages to run functions
import matplotlib
import numpy as np
import pandas
import scipy
import FlowCytometryTools
import pylab
import sys
import pathlib
import pylab
from matplotlib import pyplot as plt
from pathlib import Path
from FlowCytometryTools import test_data_dir, test_data_file
from FlowCytometryTools import FCMeasurement
from matplotlib.backends.backend_pdf import PdfPages
from FlowCytometryTools import ThresholdGate, PolyGate, QuadGate
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.patches import Rectangle
%matplotlib inline

tplate = input('What is the name of the T plate folder path?:')
sampleT = importF(tplate)
NKplate = input('What is the name of the NK plate folder path?:')
sampleNK = importF(NKplate)

treg_count = []
nonTreg_count = []
nk_count = []
bnk_count = []
cd_count = []

treg_data = []
nonTreg_data = []
nk_data = []
bnk_data = []
cd_data = []

treg_gate = treg()
nonTreg_gate = nonTreg()
nk_gate = nk()
bnk_gate = bnk()
cd_gate = cd()

for i, sample_i in enumerate(sampleT):
    treg_count.append(cellCount(sample_i, treg_gate))
    nonTreg_count.append(cellCount(sample_i, nonTreg_gate))
    
    treg_data.append(rawData(sample_i, treg_gate))
    nonTreg_data.append(rawData(sample_i, nonTreg_gate))
    
    tcells(sample_i, treg_gate, nonTreg_gate)
    
for i, sample_i in enumerate(sampleNK):
    nk_count.append(cellCount(sample_i, nk_gate))
    bnk_count.append(cellCount(sample_i, bnk_gate))
    
    nk_data.append(rawData(sample_i, nk_gate))
    bnk_data.append(rawData(sample_i, bnk_gate))
    
    nk_bnk_plot(sample_i, nk_gate, bnk_gate)
    
    
for i, sample_i in enumerate(sampleNK):
    cd_count.append(cellCount(sample_i, cd_gate))
    cd_data.append(rawData(sample_i, cd_gate))
    cd_plot(sample_i, cd_gate)
    
print('Count of T reg cells: ', treg_count)
print('Count of non T reg cells: ', nonTreg_count)

print('Count of NK cells: ', nk_count)
print('Count of Bright NK cells: ', bnk_count)
print('Count of CD cells: ', cd_count)
#print('Count of T reg cells in the file: ', treg_count)
#print('Count of non T reg cells in the file: ', nonTreg_count)
#print('Data of T reg cells: ', treg_data)
#print('Data of non T reg cells in the file: ', nonTreg_data)

# ********************** PCA Sample Program ***********************
### Import all necessary packages to run functions
import matplotlib
import numpy as np
import pandas
import scipy
import FlowCytometryTools
import pylab
import sys
import pathlib
import sklearn
from matplotlib import pyplot as plt
from pathlib import Path
from FlowCytometryTools import test_data_dir, test_data_file
from FlowCytometryTools import FCMeasurement
from matplotlib.backends.backend_pdf import PdfPages
from FlowCytometryTools import ThresholdGate, PolyGate, QuadGate
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
%matplotlib inline

import math

dataMatT = []
principalDfT = []
pstatMatT = []

dataMatNK = []
principalDfNK = []
pstatMatNK = []

from matplotlib import colors

tplate = input('What is the name of the T plate folder?:')
tsample = importF(tplate)

nkplate = input('What is the name of the NK plate folder?:')
nksample = importF(nkplate)

for i, sample_i in enumerate(tsample):
    dataT, pstatT, featuresT = sampleT(sample_i)
    dataMatT.append(dataT)
    pstatMatT.append(pstatT)
    
    xfT, loadingT = appPCA(dataT, featuresT)
    
    pcaPltCat(xfT, pstatT, loadingT, featuresT)
    
for i, sample_i in enumerate(nksample):
    dataNK, pstatNK, featuresNK = sampleNK(sample_i)
    dataMatNK.append(dataNK)
    pstatMatNK.append(pstatNK)
   
    xfNK, loadingNK = appPCA(dataNK, featuresNK)
    pcaPltCat(xfNK, pstatNK, loadingNK, featuresNK)  
   
plt.show()