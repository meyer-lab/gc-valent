"""
Theano Op for using differencing for Jacobian calculation.
"""
import numpy as np
from theano.tensor import dot, dmatrix, dvector, Op
from .model import runCkineU, runCkineSensi, nSpecies, nParams


class runCkineOp(Op):
    itypes = [dvector]
    otypes = [dvector]

    def __init__(self, ts):
        if ts.size > 1:
            raise NotImplementedError('This Op only works with a single time point.')
        
        self.ts = ts
        self.nSpecies = nSpecies()

    def infer_shape(self, node, i0_shapes):
        assert len(i0_shapes) == 1
        return [(56, )]     # TODO: this number should be different

    def perform(self, node, inputs, outputs):
        yOut, retVal = runCkineU(self.ts, inputs[0])

        assert yOut.size == self.nSpecies

        if retVal < 0:
            yOut[:] = -np.inf

        outputs[0][0] = np.squeeze(yOut)

    def grad(self, inputs, g):
        """ Calculate the runCkineOp gradient. """
        return [dot(g[0], runCkineOpDiff(self.ts)(inputs[0]))]


class runCkineOpDiff(Op):
    itypes = [dvector]
    otypes = [dmatrix]

    def __init__(self, ts):
        if ts.size > 1:
            raise NotImplementedError('This Op only works with a single time point.')

        self.ts = ts

    def perform(self, node, inputs, outputs):
        _, retVal, sensi = runCkineSensi(self.ts, inputs[0])

        if retVal < 0:
            sensi[:, :] = -np.inf

        outputs[0][0] = np.squeeze(sensi)


class runCkineKineticOp(Op):
    itypes = [dvector]
    otypes = [dvector]

    def __init__(self, ts, condense):
        self.ts = ts
        self.nSpecies = nSpecies()
        self.condense = condense
        assert condense.size == self.nSpecies

    def infer_shape(self, node, i0_shapes):
        assert len(i0_shapes) == 1
        return [(self.ts.size, )]

    def perform(self, node, inputs, outputs):
        yOut, retVal = runCkineU(self.ts, inputs[0])

        assert yOut.shape == (self.ts.size, self.nSpecies)

        if retVal < 0:
            yOut[:] = -np.inf

        outputs[0][0] = np.dot(yOut, self.condense)

    def grad(self, inputs, g):
        """ Calculate the runCkineOp gradient. """
        return [dot(g[0], runCkineOpKineticDiff(self.ts, self.condense)(inputs[0]))]


class runCkineOpKineticDiff(Op):
    itypes = [dvector]
    otypes = [dmatrix]

    def __init__(self, ts, condense):
        self.nSpecies = nSpecies()
        self.ts = ts
        self.condense = condense
        assert condense.size == self.nSpecies

    def perform(self, node, inputs, outputs):
        _, retVal, sensi = runCkineSensi(self.ts, inputs[0])

        if retVal < 0:
            sensi[:, :] = -np.inf

        outputs[0][0] = np.dot(np.transpose(sensi), self.condense)
