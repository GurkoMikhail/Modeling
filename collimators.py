from numba import jit
import numpy as np


@jit(nopython=True, cache=True)
def get_collimated(coordinates, hole_centers, hole_diameter):
    collimated = []
    for i, coord in enumerate(coordinates):
        for hole_center in hole_centers:
            if inside_hexagon(coord, hole_center, hole_diameter):
                collimated.append(i)
                break
    return collimated


@jit(nopython=True, cache=True)
def inside_hexagon(coordinates, center, d):
    x0 = coordinates[0]
    y0 = coordinates[1]
    x = center[0]
    y = center[1]
    dx = abs(x - x0)/d
    dy = abs(y - y0)/d
    a = 0.25*1.7320508
    return (dy <= a) and (a*dx + 0.25*dy <= 0.5*a)


@jit(nopython=True, cache=True)
def get_centers(size, hole_diameter, septa):
    hole_size_X = hole_diameter/2
    hole_size_Y = (hole_size_X*np.sqrt(3)/2)
    intercenter_distance_X = (hole_size_X*(1 + np.cos(np.pi/3)) + septa)*2
    intercenter_distance_Y = hole_size_Y*2 + septa
    centers = []
    x0 = hole_size_X + septa
    x1 = size[0] - (hole_size_X + septa)
    y0 = hole_size_Y + septa
    y1 = size[1] - (hole_size_Y + septa)
    i = 0
    for y in np.arange(y0, y1, intercenter_distance_Y/2):
        if i%2 == 1:
            for x in np.arange(x0, x1, intercenter_distance_X):
                centers.append((x, y))
        else:
            for x in np.arange(x0 + intercenter_distance_X/2, x1, intercenter_distance_X):
                centers.append((x, y))
        i += 1
    centers = np.asarray(centers)
    dx = (size[0] - centers[-1, 0] - (hole_size_X + septa))
    dy = (size[1] - centers[-1, 1] - (hole_size_Y + septa))
    centers[:, 1] += dy/2
    if dx > hole_size_X + septa:
        dx -= intercenter_distance_X/2
    centers[:, 0] += dx/2
    return centers


# @jit(nopython=True, cache=True)
# def get_collimated(coordinates, holes_centers, hole_diameter):
#     collimated = []
#     hole_radius2 = (hole_diameter/2)**2
#     for i, coord in enumerate(coordinates):
#         for hole_center in holes_centers:
#             r2 = (coord[0] - hole_center[0])**2 + (coord[1] - hole_center[1])**2
#             if r2 <= hole_radius2:
#                 collimated.append(i)
#                 break
#     return collimated


# def get_centers(size, hole_diameter, septa):
#     intercenter_distance = (hole_diameter + septa)
#     interlayer_distance = float((intercenter_distance/2)*np.sqrt([3]))
#     x = hole_diameter/2 + septa
#     y = hole_diameter/2 + septa
#     mindx = hole_diameter/2 + septa
#     mindy = hole_diameter/2 + septa
#     centers = []
#     layer = 0
#     len_of_layers = []
#     while y + mindy < size[1]:
#         n_in_layer = 0
#         while x + mindx < size[0]:
#             centers.append([x, y])
#             n_in_layer += 1
#             x += intercenter_distance
#         len_of_layers.append(n_in_layer)
#         y += interlayer_distance
#         if layer%2 == 0:
#             x = intercenter_distance
#         else:
#             x = intercenter_distance/2
#         layer += 1
#     max_len_of_layer = max(len_of_layers)
#     min_len_of_layer = min(len_of_layers)
#     centers = np.array(centers)
#     d = size[:2] - centers[-1] - centers[0]
#     dx = d[0]
#     if layer%2 == 0:
#         if min_len_of_layer != max_len_of_layer:
#             dx -= intercenter_distance/2
#     elif min_len_of_layer == max_len_of_layer:
#         dx -= intercenter_distance/2
#     dx /= 2
#     dy = d[1]/2
#     centers[:, 0] += dx
#     centers[:, 1] += dy
#     return centers


@jit(nopython=True, cache=True)
def inside_hexagon2(coordinates, center, d):
    x1 = coordinates[0]
    y1 = coordinates[1]
    x2 = center[0]
    y2 = center[1]
    x = np.absolute(x1 - x2)
    y = np.absolute(y1 - y2)

    py1 = d*0.86602540378
    px2 = d*0.2588190451
    py2 = d*0.96592582628

    p_angle_01 = -x*(py1 - y) - x*y
    p_angle_20 = -y*(px2 - x) + x*(py2 - y)
    p_angle_03 = y*d
    p_angle_12 = -x*(py2 - y) - (px2 - x)*(py1 - y)
    p_angle_32 = (d - x)*(py2 - y) + y*(px2 - x)

    is_inside_1 = (p_angle_01*p_angle_12 >= 0) and (p_angle_12*p_angle_20 >= 0)
    is_inside_2 = (p_angle_03*p_angle_32 >= 0) and (p_angle_32*p_angle_20 >= 0)

    return is_inside_1 or is_inside_2
    # return np.nonzero(is_inside_1 + is_inside_2)[0]


if __name__ == "__main__":

    # from pyqtgraph.Qt import QtGui, QtCore
    import numpy as np
    # import pyqtgraph as pg

    centers = get_centers((2., 2.), 0.15, 0.02)
    print(centers.shape[0])

    # app = QtGui.QApplication([])
    # win = QtGui.QMainWindow()
    # win.resize(1920,1080)
    # imv = pg.ImageView()
    # win.setCentralWidget(imv)
    # win.show()
    # win.setWindowTitle('ImageCollimator')

    coordinates = []
    for x in np.linspace(0, 2., 10**3):
        for y in np.linspace(0, 2., 10**3):
            coordinates.append((x, y))
    coordinates = np.array(coordinates)
    collimated = get_collimated(coordinates, centers, 0.15)
    coordinates = np.delete(coordinates, collimated, axis=0)
    x = coordinates[:, 0]
    y = coordinates[:, 1]

    histogram = np.histogram2d(x, y, bins=int(10**3))
    data = histogram[0]

    import pyqtgraph as pg
    from pyqtgraph.Qt import QtCore, QtGui


    # Interpret image data as row-major instead of col-major
    pg.setConfigOptions(imageAxisOrder='row-major')

    pg.mkQApp()
    win = pg.GraphicsLayoutWidget()
    win.setWindowTitle('pyqtgraph example: Image Analysis')

    # A plot area (ViewBox + axes) for displaying the image
    p1 = win.addPlot(title="")

    # Item for displaying image data
    img = pg.ImageItem()
    p1.addItem(img)

    # Custom ROI for selecting an image region
    roi = pg.ROI([0, 14], [6, 5])
    roi.addScaleHandle([0.5, 1], [0.5, 0.5])
    roi.addScaleHandle([0, 0.5], [0.5, 0.5])
    p1.addItem(roi)
    roi.setZValue(10)  # make sure ROI is drawn above image

    # Isocurve drawing
    iso = pg.IsocurveItem(level=0.8, pen='g')
    iso.setParentItem(img)
    iso.setZValue(5)

    # Contrast/color control
    hist = pg.HistogramLUTItem()
    hist.setImageItem(img)
    win.addItem(hist)

    # Draggable line for setting isocurve level
    isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
    hist.vb.addItem(isoLine)
    hist.vb.setMouseEnabled(y=False) # makes user interaction a little easier
    isoLine.setValue(0.8)
    isoLine.setZValue(1000) # bring iso line above contrast controls

    # Another plot area for displaying ROI data
    win.nextRow()
    p2 = win.addPlot(colspan=2)
    p2.setMaximumHeight(250)
    win.resize(800, 800)
    win.show()


    # Generate image data
    # data = np.random.normal(size=(200, 100))
    # data[20:80, 20:80] += 2.
    # data = pg.gaussianFilter(data, (3, 3))
    # data += np.random.normal(size=(200, 100)) * 0.1
    img.setImage(data)
    hist.setLevels(data.min(), data.max())

    # build isocurves from smoothed data
    iso.setData(pg.gaussianFilter(data, (2, 2)))

    # set position and scale of image
    img.scale(0.02, 0.02)
    # img.translate(-50, 0)

    # zoom to fit imageo
    p1.autoRange()  


    # Callbacks for handling user interaction
    def updatePlot():
        global img, roi, data, p2
        selected = roi.getArrayRegion(data, img)
        p2.plot(selected.mean(axis=0), clear=True)

    roi.sigRegionChanged.connect(updatePlot)
    updatePlot()

    def updateIsocurve():
        global isoLine, iso
        iso.setLevel(isoLine.value())

    isoLine.sigDragged.connect(updateIsocurve)

    def imageHoverEvent(event):
        """Show the position, pixel, and value under the mouse cursor.
        """
        if event.isExit():
            p1.setTitle("")
            return
        pos = event.pos()
        i, j = pos.y(), pos.x()
        i = int(np.clip(i, 0, data.shape[0] - 1))
        j = int(np.clip(j, 0, data.shape[1] - 1))
        val = data[i, j]
        ppos = img.mapToParent(pos)
        x, y = ppos.x(), ppos.y()
        p1.setTitle("pos: (%0.1f, %0.1f)  pixel: (%d, %d)  value: %g" % (x, y, i, j, val))

    # Monkey-patch the image to use our custom hover function. 
    # This is generally discouraged (you should subclass ImageItem instead),
    # but it works for a very simple use like this. 
    img.hoverEvent = imageHoverEvent
    QtGui.QApplication.instance().exec_()
