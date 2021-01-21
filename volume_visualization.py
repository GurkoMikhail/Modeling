from PyQt5 import QtGui
from PyQt5.uic import loadUi
from pyqtgraph.Qt import mkQApp
from pyqtgraph.graphicsItems.HistogramLUTItem import HistogramLUTItem
from pyqtgraph.opengl.items.GLAxisItem import GLAxisItem
from pyqtgraph.opengl.items.GLBoxItem import GLBoxItem

mkQApp()
mainWindow = loadUi("UI/volume_visualization.ui")

hist = HistogramLUTItem()
# hist.setImageItem(image)
mainWindow.graphicsView_histogram.addItem(hist)

space_box = GLBoxItem()
space_box.setSize(10, 10, 10)
mainWindow.graphicsView_volume.addItem(space_box)


mainWindow.show()
QtGui.QApplication.exec()