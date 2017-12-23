import theano.tensor as T
import numpy as np


class centralDiff(T.Op):
    itypes = [T.dvector]
    otypes = [T.dvector]

    def __init__(self, calcModel):
        self.M = calcModel
        self.dg = centralDiffGrad(calcModel)

    def infer_shape(self, node, i0_shapes):
        return [(self.M.concs, )]

    def perform(self, node, inputs, outputs):
        vec, = inputs
        mu = self.M.calc(vec)

        outputs[0][0] = np.array(mu)

   # def grad(self, inputs, g):
   #     return T.dot(g, self.dg(inputs[0]))


class centralDiffGrad(T.Op):
    itypes = [T.dvector]
    otypes = [T.dmatrix]

    def __init__(self, calcModel):
        self.M = calcModel

    def infer_shape(self, node, i0_shapes):
        return [(self.M.concs, )]

    def num_grad(self, inputs, g):
        # Find our current point
        x0 = inputs[0]
        f0 = self(x0)
        epsilon = 1.0E-7

        print(T.shape(x0))

        jac = T.dot(x0, f0)
        dx = T.zeros_like(x0)

        for i in range(len(x0)):
            dx[i] = epsilon
            jac[i] = (self(x0 + dx) - f0)/epsilon
            dx[i] = 0.0

        return T.dot(g, jac.T)
    