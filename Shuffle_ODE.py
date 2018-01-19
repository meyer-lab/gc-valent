from ckine.model import solveAutocrine, fullModel, __active_species_IDX
import numpy as np


def approx_jacobian():
    """Approximate the Jacobian matrix of callable function func

       * Parameters
         x       - The state vector at which the Jacobian matrix is desired
         func    - A vector-valued function of the form f(x,*args)
         epsilon - The peturbation used to determine the partial derivatives
         *args   - Additional arguments passed to func

       * Returns
         An array of dimensions (lenf, lenx) where lenf is the length
         of the outputs of func, and lenx is the number of

       * Notes
         The approximation is done using forward differences

    """
    active_species_IDX = __active_species_IDX
    trafRates = np.random.sample(17)
    #t = np.array([0.0, 100000.0])
    rxnRates = np.random.sample(17)
    #args = (t, rxnRates, trafRates, __active_species_IDX)
    x0 = np.random.sample(56)
    f0 = fullModel(x0, 0.0 , rxnRates, trafRates, active_species_IDX)
    jac = np.zeros([len(x0),len(f0)])
    dx = np.zeros(len(x0))
    for i in range(len(x0)):
       dx[i] = np.sqrt(np.finfo(float).eps)
       jac[i] = (fullModel(x0 + dx, 0.0 , rxnRates, trafRates, active_species_IDX) - f0)/(np.sqrt(np.finfo(float).eps))
       dx[i] = 0.0
    return jac.transpose()

a = approx_jacobian()
np.savetxt('Nonzero Boolean', a, delimiter=',')
print(a)
print(a != 0)