from ckine.model import solveAutocrine, fullModel, __active_species_IDX
import numpy as np
import scipy


def approx_jacobian():
    """Approximate the Jacobian matrix of callable function func

       * Parameters
         None

       * Returns
         An array of dimensions (lenf, lenx) where lenf is the length
         of the outputs of func, and lenx is the number of

       * Notes
         The approximation is done using forward differences. func in this case is the fullModel function from the main model file. 

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
np.savetxt('Nonzero Boolean.csv', (a != 0).astype(np.int), fmt='%d', delimiter=' ')



#b = a!=0

#bm= scipy.sparse.csr_matrix(b)
#permb = scipy.sparse.csgraph.reverse_cuthill_mckee(bm, False)
#B = bm[np.ix_(permb,permb)].A
#np.savetxt('Shuffled Nonzero Boolean.csv', B, fmt='%d', delimiter=' ')
