from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget
from PyQt5.uic import loadUi
from numpy import array, asarray, load, min, max, arange, nonzero, uint32, float32, NaN
from pyqtgraph.Qt import mkQApp
from pyqtgraph import GradientEditorItem, makeRGBA
from pyqtgraph.opengl import GLBoxItem, GLVolumeItem

voxel_size = 0.4

volume = load('Phantoms/efg3cut_source.npy').astype(float32)
volume_size = array(volume.shape, dtype=float32)*voxel_size

levels = [
    min(volume),
    max(volume)
]

def gradientChanged():
    global volume, volumeItem, gradientEditor, levels
    listTicks = gradientEditor.listTicks()
    listTicks[0][0].color.setAlpha(0)
    for tick in listTicks[1:]:
        tick[0].color.setAlpha(30)
    lut = gradientEditor.getLookupTable(255)
    volume_colors = asarray([makeRGBA(data=slice, lut=lut, levels=levels)[0] for slice in volume])
    volumeItem.setData(volume_colors)

mkQApp()
mainWindow = loadUi("UI/volume_visualization.ui")

volumeItem = GLVolumeItem(None, sliceDensity=2)
volumeItem.scale(*[voxel_size]*3)
gradientEditor = GradientEditorItem(orientation='right')
gradientEditor.sigGradientChangeFinished.connect(gradientChanged)
mainWindow.graphicsView.addItem(gradientEditor)
volumeViewWidget = mainWindow.openGLWidget
volumeViewWidget.addItem(volumeItem)

volumeViewWidget.setCameraPosition(distance=volume_size[0]*3)
volumeViewWidget.pan(*volume_size/2)

space_box = GLBoxItem()
space_box.setSize(*volume_size)
mainWindow.openGLWidget.addItem(space_box)


mainWindow.show()
QtGui.QApplication.exec()