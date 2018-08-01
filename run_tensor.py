"""File to pickle factorization into a new file."""
import os
import pickle
from ckine.tensor_generation import prepare_tensor
from ckine.Tensor_analysis import perform_decomposition

values, _, _, _, cell_names = prepare_tensor(2)
factors = perform_decomposition(values, 20)

filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), "./ckine/data/factors_results/Sampling.pickle")
with open(filename, 'wb') as f:
    pickle.dump((factors, cell_names), f)