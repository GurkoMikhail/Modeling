from numba import vectorize, float64
from math import log, exp, nan, sqrt, acos, pi, cos
from random import random

from numpy.core.numeric import NaN

@vectorize([float64(float64)], nopython=True, cache=True)
def generation_theta(energy):
    E0_m = energy/510998.9461
    eps0 = 1/(1 + 2*E0_m)
    epsilon0sq = eps0*eps0
    alpha1 = -log(eps0)
    alpha2 = (1 - epsilon0sq)/2
    onecost = 0.
    greject = 0.
    while greject < random():
        if alpha1/(alpha1 + alpha2) > random():
            epsilon = exp(-alpha1*random())
            epsilonsq = epsilon*epsilon
        else:
            epsilonsq = epsilon0sq + (1 - epsilon0sq)*random()
            epsilon = sqrt(epsilonsq)
        onecost = (1 - epsilon)/(epsilon*E0_m)
        sint2 = onecost*(2 - onecost)
        greject = 1 - epsilon*sint2/(1+ epsilonsq)
    costheta = 1 - onecost
    theta = acos(costheta)
    return theta

@vectorize([float64(float64, float64, float64, float64)], nopython=True, cache=True)
def culculate_sigma(energy, Z, density, M):
    mec2 = 510998.9461          #eV
    re = 2.8179403267*10**(-13) #cm
    Na = 6.02214076*10**23      #mol^(-1)
    k = energy/mec2
    k_1 = 1 + k
    k2 = 2*k
    k2_1 = 1 + k2
    ln1_2k = log(k2_1)
    sigma = 2*pi*(re**2)*Z
    sigma *= (k_1/(k**2))*(2*k_1/k2_1 - ln1_2k/k) + ln1_2k/k2 - (1 + 3*k)/(k2_1**2)
    sigma *= density
    sigma *= Na
    sigma /= M
    return sigma

@vectorize([float64(float64, float64)], nopython=True, cache=True)
def culculate_energy_change(energy, theta):
    mec2 = 510998.9461          #eV
    k = energy/mec2
    k1_cos = k*(1 - cos(theta))
    energy_change = energy*k1_cos/(1+k1_cos)
    return energy_change

if __name__ == "__main__":
    # import numpy as np
    # from time import time
    # energy = np.linspace(100*10**3, 140*10**3, 10**3, dtype=np.float64)
    # print ('Start!')
    # start = time()
    # mu = culculate_sigma(energy, 10, 1, 18)
    # finish = time() - start
    # print(f'Finish!\n{finish} seconds')

    # from pyqtgraph.Qt import QtCore, QtGui
    # import pyqtgraph as pg
    # app = QtGui.QApplication([])
    # pg.setConfigOption('background', 'w')
    # pg.setConfigOption('foreground', 'k')
    # mw = QtGui.QMainWindow()
    # mw.resize(1920,1080)
    # view = pg.GraphicsLayoutWidget()
    # mw.setCentralWidget(view)
    # mw.setWindowTitle('Compton')
    # plt = view.addPlot()
    # scatter = pg.ScatterPlotItem(x=energy/10**6, y=mu, size=2)
    # plt.addItem(scatter)
    # mw.show()
    # QtGui.QApplication.exec()

    import numpy as np
    from time import time
    energy = np.full(10**8, 140*10**6, dtype=np.float64)
    indices = np.array([0, 1, 2, 3, 4, 5, np.inf, 7])
    # energy = energy[indices]
    # energy = np.array([np.inf], dtype=np.float64)
    # energy = 7
    print ('Start!')
    start = time()
    theta = generation_theta(energy)
    finish = time() - start
    print(f'Finish!\n{finish} seconds')

    hist = np.histogram(np.cos(theta), 10**3, density=False)

    from pyqtgraph.Qt import QtCore, QtGui
    import pyqtgraph as pg
    app = QtGui.QApplication([])
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    mw = QtGui.QMainWindow()
    mw.resize(1920,1080)
    view = pg.GraphicsLayoutWidget()
    mw.setCentralWidget(view)
    mw.setWindowTitle('Compton')
    plt = view.addPlot()
    scatter = pg.ScatterPlotItem(x=np.arccos(hist[1][:-1]), y=hist[0], size=2)
    plt.addItem(scatter)
    mw.show()
    QtGui.QApplication.exec()
