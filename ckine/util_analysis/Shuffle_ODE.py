from ckine.model import fullModel, dy_dt
import numpy as np


def approx_jacobian(y, t, rxn, tfr, delta=1.0E-9):
    """Approximate the Jacobian matrix of callable function func

       * Parameters
         None

       * Returns
         An array of dimensions (lenf, lenx) where lenf is the length
         of the outputs of func, and lenx is the number of

       * Notes
         The approximation is done using forward differences. func in this case is the fullModel function from the main model file. 

    """
    f0 = fullModel(y, t, rxn, tfr)
    jac = np.zeros([len(y),len(f0)])
    dy = np.zeros(len(y))
    for i in range(y.size):
        dy[i] = delta
        jac[i] = (fullModel(y + dy, t, rxn, tfr) - f0)/delta
        dy[i] = 0.0
        
    return jac.transpose()


def approx_jac_dydt(y, t, rxn, delta=1.0E-9):
    """Approximate the Jacobian matrix of callable function func

       * Returns
         An array of dimensions (lenf, lenx) where lenf is the length
         of the outputs of func, and lenx is the number of

       * Notes
         The approximation is done using forward differences. func in this case is the fullModel function from the main model file. 

    """
    f0 = dy_dt(y, t, rxn)
    jac = np.zeros([y.size, f0.size])

    for i in range(y.size):
        dy = y.copy()
        dy[i] = dy[i] + delta

        jac[i] = (dy_dt(dy, t, rxn) - f0)/delta
        
    return jac.transpose()
