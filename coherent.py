import numpy as np
from materials import Materials
from numpy import arccos, arcsin, exp, cos, sin, pi
from h5py import File
from scipy.interpolate import interp2d 

re = 2.8179403267*10**(-15) #m
re2 = (re**2)*10**28        #barn
hd = 6.582199514*10**(-16)*2*pi
cLight = 299792458

class Indicatrix:

    def __init__(self, meterials):
        self.meterials = meterials
        self.base = 'PHY.F20'
        self.constructTable()

    def constructTable(self):
        self.table = {}
        file = File('tables/form-factors.hdf')
        baseGroup = file[f'{self.base}']
        for index, material in self.meterials.materials_dict.items():
            coefficientsGroup = baseGroup[material]
            a = np.copy(coefficientsGroup['a'])
            b = np.copy(coefficientsGroup['b'])
            c = np.copy(coefficientsGroup['c'])
            self.table.update({
                index: (a, b, c)
            })
        file.close()

    def atomicFormFactor(self, material, q):
        a, b, c = self.table[material]
        formFactor = np.full_like(q, c)
        for i in range(4):
            formFactor += a[i]*exp(-b[i]*(q/(4*pi))**2)
        return formFactor

    def value_(self, material, theta, energy):
        # lamb = (10**10)*hd*cLight/energy
        q = (sin(theta/2)*energy/12398.520)
        # q = 2*energy*sin(theta/2)/30000
        # q = 2*sin(theta/2)*energy/12398.570430885673
        p = re2/2*(1 + cos(theta)**2)*(self.atomicFormFactor(material, q))**2
        return p
        
    def value(self, material, x, energy):
        # lamb = (10**10)*hd*cLight/energy
        theta = arccos(x)
        q = sin(theta/2)*energy/12398.520
        # q = 2*energy*sin(theta/2)/30000
        # q = 2*sin(theta/2)*energy/12398.570430885673
        p = re2/2*(1 + x**2)*(self.atomicFormFactor(material, q))**2
        return p



if __name__ == '__main__':
    materials = {
        'Elemental media/Pb':                                       0,
    }
    materials = Materials(materials, max_energy=140500.)
    indicatrix = Indicatrix(materials)
    cosTheta = np.linspace(-1, 1, 50)
    theta = np.linspace(0, pi, 50)
    energy = np.linspace(1000, 100000, 50)
    cosTheta_, energy_ = np.stack(np.meshgrid(cosTheta, energy), axis=0).reshape((2, -1))
    value = indicatrix.value(0, cosTheta_, energy_)

    func = interp2d(value, energy_, cosTheta_)
    # func = interp2d(arccos(cosTheta_), energy_, value)

    import pyqtgraph as pg
    from PyQt5 import QtGui

    pg.mkQApp()
    # pg.plot(cos(theta))

    # ksi = np.linspace(value.min(), value.max(), 50)

    # data = func(theta, energy)
    data = arccos(func(value, energy))
    # data = value.reshape((energy.size, theta.size))

    imv = pg.ImageView()
    imv.setImage(data)
    imv.show()
    QtGui.QApplication.instance().exec_()


