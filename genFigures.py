#!/usr/bin/env python3

from ckine.figures.figureCommon import overlayCartoon
import sys
import logging
import time
import matplotlib
matplotlib.use('AGG')

fdir = './output/'

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

if __name__ == '__main__':
    start = time.time()
    nameOut = 'figure' + sys.argv[1]

    exec('from ckine.figures import ' + nameOut)
    ff = eval(nameOut + '.makeFigure()')
    ff.savefig(fdir + nameOut + '.svg', dpi=ff.dpi, bbox_inches='tight', pad_inches=0)

    if sys.argv[1] == 'C2':
        # Overlay Figure 2 cartoon
        overlayCartoon(fdir + 'figureC2.svg',
                       './ckine/graphics/tensor4D.svg', 5, 6, scalee=1.62)

    if sys.argv[1] == 'C3':
        # Overlay Figure 3 cartoon
        overlayCartoon(fdir + 'figureC3.svg',
                       './ckine/graphics/ModelCartoon.svg', 1450, 0, scalee=0.023)

    logging.info('%s is done after %s seconds.', nameOut, time.time() - start)
