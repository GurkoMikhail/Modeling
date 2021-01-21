import h5py
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
from scipy.ndimage import gaussian_filter
import pyqtgraph as pg
import time as tm

time = 30.
detector_size = np.array([51.2, 40, 3.])
pixel_size = 0.1
resolution = 0.4/2.355

range_size = np.column_stack([[0, 0], detector_size[:2]])

coordinates = []
energy_transfer = []
emission_time = []
emission_coordinates = []

start = tm.time()

name = '<modeling.Modeling object at 0x7fbfe8304df0>'

file = h5py.File(f'Output data/{name}', 'r')
for group in file['Flows']:
    flow = file['Flows'][group]
    coordinates.append(flow['Coordinates'])
    energy_transfer.append(flow['Energy transfer'])
    emission_time.append(flow['Emission time'])
    emission_coordinates.append(flow['Emission coordinates'])
    # time = group[-6:-1]
coordinates = np.concatenate(coordinates)
energy_transfer = np.concatenate(energy_transfer)
emission_time = np.concatenate(emission_time)
emission_coordinates = np.concatenate(emission_coordinates)
file.close()

print(tm.time() - start)

def indices_repeated(emission_time):
    unique, unique_counts = np.unique(emission_time, return_counts=True)
    indices = np.nonzero(unique_counts > 1)[0]
    unique = unique[indices]
    indices = []
    for time in unique:
        indices.append(np.nonzero(emission_time == time)[0])
    return indices

def average_coordinates(coordinates, energy, emission_time):
    indices = indices_repeated(emission_time)
    del_indices = []
    for event_indices in indices:
        weights = energy[event_indices]
        energy[event_indices[0]] = np.sum(weights)
        average = np.average(coordinates[event_indices], axis=0, weights=weights)
        coordinates[event_indices[0]] = average
        del_indices.append(event_indices[1:])
    del_indices = np.concatenate(del_indices)
    coordinates = np.delete(coordinates, del_indices)
    energy = np.delete(energy, del_indices)
    emission_time = np.delete(emission_time, del_indices)

def energy_deviation(energy, coeff):
    resolution_distribution = coeff/np.sqrt(energy/10**6)
    sigma = resolution_distribution*energy_transfer/2.355
    energy[:] = np.random.normal(energy, sigma)

def in_subject(min_level, max_level):
    return (coordinates[:, 1] >= min_level)*(coordinates[:, 1] <= max_level)

def enough_energy(min_energy, max_energy):
    return (energy_transfer >= min_energy)*(energy_transfer <= max_energy)

def updateImage():
    global image, data, p2
    min_energy = minEnergyLine.value()
    max_energy = maxEnergyLine.value()
    min_level = 0
    max_level = 3
    indices = np.nonzero(enough_energy(min_energy, max_energy))[0]
    histogram = np.histogram2d(coordinates[indices, 0], coordinates[indices, 2], bins=(detector_size[:2]/pixel_size).astype(np.int)[:2], range=range_size[:2])
    data = histogram[0][::-1]
    # data = np.rot90(data)
    data = gaussian_filter(data, (resolution/pixel_size, resolution/pixel_size))
    image.setImage(data)
    hist.setLevels(data.min(), data.max())

inside = in_subject(0., 3.)
coordinates = coordinates[inside]
energy_transfer = energy_transfer[inside]
emission_time = emission_time[inside]
average_coordinates(coordinates, energy_transfer, emission_time)
coeff = np.sqrt(0.14)*0.099
energy_deviation(energy_transfer, coeff)


pg.mkQApp()
win = pg.GraphicsLayoutWidget()
win.setWindowTitle(f'Image in {time} seconds')

p1 = win.addPlot()
p1.setLabel('left', "Y Axis", units='sm')
p1.setLabel('bottom', "X Axis", units='sm')

image = pg.ImageItem()
p1.addItem(image)
image.scale(pixel_size, pixel_size)

hist = pg.HistogramLUTItem()
hist.setImageItem(image)
win.addItem(hist)

win.nextRow()
p2 = win.addPlot(colspan=2, title='Energy distribution')
p2.setLabel('left', 'N')
p2.setLabel('bottom', 'Energy', units='eV')
p2.setLogMode(y=True)
p2.setMaximumHeight(250)
win.resize(1200, 1000)
win.show()
energy_distribution, energy = np.histogram(energy_transfer, bins=256)
energy_distribution += 1

p2.plot(x=energy[1:], y=energy_distribution, clear=True)

minEnergyLine = pg.InfiniteLine(120*10**3, movable=True, markers=[('^', 0), ('v', 1)], pen='y')
p2.addItem(minEnergyLine)
minEnergyLine.sigPositionChangeFinished.connect(updateImage)

maxEnergyLine = pg.InfiniteLine(160*10**3, movable=True, markers=[('^', 0), ('v', 1)], pen='y')
p2.addItem(maxEnergyLine)
maxEnergyLine.sigPositionChangeFinished.connect(updateImage)

updateImage()

QtGui.QApplication.instance().exec_()