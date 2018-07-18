"""
This creates Figure S2.
"""
import string
from .figureCommon import subplotLabel, getSetup
from ..tensor_generation import prepare_tensor
from ..Tensor_analysis import perform_decomposition
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    x, y = 3, 5
    ax, f = getSetup((10, 10), (x, y))

    factors, new_mat, mat, mats, cell_names = decompose_tensor(2)
    a = b = c = d = -1
    for ii in range(x*y):
        if ii % y == 0:
            a += 2
            plot_timepoint(ax[ii],factors,a, a+1)
            
        elif ii % y == 1:
            b += 2
            plot_cells(ax[ii],factors, mat, b, b+1, cell_names)
            
        elif ii % y == 2:
            c+=2
            plot_ligands(ax[ii],factors,c, c+1)

        elif ii % y == 3:
            d+= 2
            plot_values(ax[ii],factors,d, d+1)

    f.tight_layout()

    return f

def decompose_tensor(lig, n_components = 6):
    """Generate the tensor and decompose it using tensor factorization."""
    #lig is how many values of ligand stimulation to use. n_components is how many components to decompose the tensor into. 
    values, new_mat, mat, mats, cell_names = prepare_tensor(lig)
    factors = perform_decomposition(values, n_components)
    return factors, new_mat, mat, mats, cell_names

def plot_values(ax,factors,component_x, component_y):
    """Plot the values decomposition factors matrix."""
    #Generate a plot for component x vs component y of the factors[3] above representing our values
    labels = ['IL2', 'IL15', 'IL7', 'IL9', 'IL4','IL21','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra','IL2Ra', 'IL2Rb', 'gc', 'IL15Ra', 'IL7Ra', 'IL9R', 'IL4Ra','IL21Ra']
    #Set Active to color red. Set Surface to color blue. Set Total to color black
    for i in range(len(factors[3])):
        data = (factors[3][:, component_x - 1][i], factors[3][:, component_y - 1][i])
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

    set_axis_limits(ax, factors, 3, component_x, component_y)

def set_axis_limits(ax, factors, f_pos, component_x, component_y):
    """Function to set axes labels and limits of the designated factor matrix."""
    ax.set_xlabel('Component ' + str(component_x))
    ax.set_ylabel('Component ' + str(component_y))
    ax.set_xlim(-np.max(np.absolute(factors[f_pos][:,component_x - 1]))*1.05,np.max(np.absolute(factors[f_pos][:,component_x - 1]))*1.05)
    ax.set_ylim(-np.max(np.absolute(factors[f_pos][:,component_y - 1]))*1.05,np.max(np.absolute(factors[f_pos][:,component_y - 1]))*1.05)
    
def plot_timepoint(ax,factors,component_x, component_y):
    """Plot the timepoint decomposition in the first column of figS2."""
    for i in range(len(factors[0])):
        ax.scatter(factors[0][:,component_x - 1][i], factors[0][:,component_y - 1][i], color = 'k')
        if i == 999:
            ax.annotate(str(i+1), xy=(factors[0][:,component_x - 1][i], factors[0][:,component_y - 1][i]), xytext = (0, 0), textcoords = 'offset points')
    set_axis_limits(ax, factors, 0, component_x, component_y)

def plot_cells(ax,factors, mat, component_x, component_y, cell_names):
    """This function plots the combination decomposition based on cell type."""
    colors = cm.rainbow(np.linspace(0, 1, len(cell_names)))
    ax.scatter(factors[1][:,component_x - 1], factors[1][:,component_y - 1], c=colors, label = cell_names)
    set_axis_limits(ax, factors, 1, component_x, component_y)

def plot_ligands(ax, factors, component_x, component_y):
    "This function is to plot the ligand combination dimension of the values tensor."
    ax.scatter(factors[2][:,component_x - 1], factors[2][:,component_y - 1])
    set_axis_limits(ax, factors, 2, component_x, component_y)
