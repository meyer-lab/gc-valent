"""
Analyze tensor from tensor_generation and plotting.
"""
import numpy as np
import tensorly as tl
from tensorly.decomposition import parafac, tucker, non_negative_tucker
from tensorly.metrics.regression import variance as tl_var, standard_deviation

backend = 'numpy' # Make backend setting accessible if needed
tl.set_backend(backend) # Set the backend within every file that imports this file


def R2X(reconstructed, original):
    ''' Calculates R2X of two tensors. '''
    return 1.0 - tl_var(reconstructed - original) / tl_var(original)

def perform_parafac(tensor, r):
    '''Apply z scoring and perform PARAFAC decomposition'''
    return parafac(tensor, rank=r) # can do verbose and tolerance (tol)

def perform_tucker(tensor, rank_list):
    '''Function to peform tucker decomposition.'''
    return non_negative_tucker(tensor, rank=rank_list, n_iter_max=1000, tol=1.0E-6) # out is (core tensor, factors)

def find_R2X_tucker(values, out):
    '''Compute R2X for the tucker decomposition.'''
    return R2X(tl.tucker_to_tensor(out[0], out[1]), values)

def find_R2X_parafac(values, factors):
    '''Compute R2X. Note that the inputs values and factors are in numpy.'''
    return R2X(tl.kruskal_to_tensor(factors), values)

def reorient_one(factors, component_index):
    """Function that takes in the 4 factor matrices and decides if that column index should flip or not and then flips it."""
    factors_idx = [factors[0][:,component_index], factors[1][:,component_index], factors[2][:,component_index]]
    component_means = tl.tensor([tl.mean(factors_idx[0]**3), tl.mean(factors_idx[1]**3), tl.mean(factors_idx[2]**3)])
    if tl.sum(component_means < 0) >= 2 and tl.sum(component_means < 0) <= 3: #if at least 2 are negative, then flip one
        count = 1
        for index, factor_idx in enumerate(factors_idx):
            if component_means[index] < 0 and count < 3:
                factors[index][:, component_index] = factor_idx * -1
                count += 1
    return factors

def reorient_factors(factors):
    """This function is to reorient the factors if at least one component in two factors matrices are negative."""
    for jj in range(factors[0].shape[1]):
        factors = reorient_one(factors, jj)
    return factors

def scale_time_factors(factors, component_index):
    """Scale the timepoint factor component by dividing the mean and then in the values plot multiply the values by that same number."""
    scale_factor = tl.mean(factors[0][:, component_index])
    factors[2][:, component_index] *= scale_factor
    factors[0][:, component_index] /= scale_factor
    return factors

def scale_all(factors):
    """Function to rescale all components. Timepoint factor matrix and values factor matrix."""
    for ii in range(factors[0].shape[1]):
        factors = scale_time_factors(factors, ii)
    return factors
