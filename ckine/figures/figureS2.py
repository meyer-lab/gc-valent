"""
This creates Figure S2.
"""
import string
from .figureCommon import subplotLabel, getSetup
from ..tensor_generation import findy, reduce_values
from ..Tensor_analysis import perform_decomposition
import numpy as np
import matplotlib.pyplot as plt

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((7, 6), (3, 4))

    factors = decompose_tensor(2)
    counter = 0
    for ii in range(3*4):
        if ii % 4 == 0:
            counter += 1
            plot_values(ax[ii],factors,counter, counter+1)
    
    # Add subplot labels
    #for ii, item in enumerate(ax):
        #subplotLabel(item, string.ascii_uppercase[ii])

    f.tight_layout()

    return f

def decompose_tensor(lig, n_components = 6):
    """Generate the tensor and decompose it using tensor factorization."""
    #lig is how many values of ligand stimulation to use. n_components is how many components to decompose the tensor into. 
    y_combos = findy(lig)[0]
    values = reduce_values(y_combos)
    factors = perform_decomposition(values, n_components)
    return factors

def plot_values(ax,factors,component_x, component_y):
    """Plot the values decomposition factors matrix."""
    #Generate a plot for component x vs component y of the factors[2] above representing our values
    labels = ['IL2', 'IL15', 'IL7', 'IL9', 'IL4','IL21','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra']
    #Set Active to color red. Set Surface to color blue. Set Total to color black
    for i in range(len(factors[2])):
        data = (factors[2][:, component_x - 1][i], factors[2][:, component_y - 1][i])
        if i in range(6):
            c = 'r'
            if i==0:
                ax.scatter(data[0], data[1], color = c, label = 'Ligand Activity')
            else:
                ax.scatter(data[0], data[1], color = c)
            ax.annotate(labels[i], xy=data, xytext = (0, 0), textcoords = 'offset points')
        elif i in range(6, 14):
            c = 'b'
            if i == 6:
                ax.scatter(data[0], data[1], color = c, label = 'Surface Receptor')
            else:
                ax.scatter(data[0], data[1], color = c)
            ax.annotate(labels[i], xy=data, xytext = (0, 0), textcoords = 'offset points')
        else:
            c = 'k'
            if i==14:
                ax.scatter(data[0], data[1], color = c, label = 'Total Receptor')
            else:
                ax.scatter(data[0], data[1], color = c)
            ax.annotate(labels[i], xy=data, xytext = (0, 0), textcoords = 'offset points')

    ax.set_xlabel('Component ' + str(component_x))
    ax.set_ylabel('Component ' + str(component_y))
    ax.set_title('Values decomposition')
    #ax.legend()    
    
