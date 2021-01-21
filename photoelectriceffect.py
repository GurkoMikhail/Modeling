from numba import vectorize, float64
from math import cos, pi, log, sqrt

@vectorize([float64(float64, float64, float64, float64)], nopython=True, cache=True)
def culculate_sigma(energy, Z, density, M):
    mec2 = 510998.9461          #eV
    re = 2.8179403267*10**(-13) #cm
    Na = 6.02214076*10**23      #mol^(-1)
    alpha = 1/137.03598
    k = energy/mec2
    sigma = 4*pi*re**2
    sigma *= alpha**4
    sigma *= Z**5
    sigma *= (2 + k)**(3/2)
    sigma /= k**(7/2)
    sigma *= 4/3 + ((1 + k)*(k - 1)/(2 + k))*(1 - log((1 + k + sqrt(k*(2 + k)))/(1 + k - sqrt(k*(2 + k)))))
    sigma *= density
    sigma *= Na
    sigma /= M
    return sigma


if __name__ == "__main__":
    import numpy as np
    from time import time
    energy = np.linspace(20, 140, 10**3)
    print ('Start!')
    start = time()
    sigma = culculate_sigma(energy, 7, 1.2*10**(-6), 14)
    finish = time() - start
    print(f'Finish!\n{finish} seconds')

    # hist = np.histogram(np.cos(theta), 10**3, density=False)

    from pyqtgraph.Qt import QtCore, QtGui
    import pyqtgraph as pg
    app = QtGui.QApplication([])
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')
    mw = QtGui.QMainWindow()
    mw.resize(1920,1080)
    view = pg.GraphicsLayoutWidget()
    mw.setCentralWidget(view)
    mw.setWindowTitle('PhotoelectircEffect')
    plt = view.addPlot()
    scatter = pg.ScatterPlotItem(x=energy, y=sigma, size=2)
    plt.addItem(scatter)
    mw.show()
    QtGui.QApplication.exec()