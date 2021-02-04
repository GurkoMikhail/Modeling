from numba import jit
import numpy as np
from math import sqrt, cos, sin, pi
from random import random

from numpy.lib.function_base import vectorize


# @jit(nopython=True, cache=True)
def culculate_R(euler_angles):
    alpha, beta, gamma = euler_angles

    Rx = np.array([
        [1,                 0,              0               ],
        [0,                 np.cos(beta),   np.sin(beta)    ],
        [0,                 -np.sin(beta),  np.cos(beta)    ]
    ])

    Ry = np.array([
        [np.cos(gamma),     0,              -np.sin(gamma)  ],
        [0,                 1,              0               ],
        [np.sin(gamma),     0,              np.cos(gamma)   ]
    ])

    Rz = np.array([
        [np.cos(alpha),     np.sin(alpha),  0               ],
        [-np.sin(alpha),    np.cos(alpha),  0               ],
        [0,                 0,              1               ]
    ])

    Rzx = np.dot(Rz, Rx)
    R = np.dot(Rzx, Ry)
    return R


@jit(nopython=True, cache=True)
def rotating_the_coordinates(coordinates, R):
    for i in range(coordinates.shape[0]):
        np.dot(R, coordinates[i], out=coordinates[i])


# @jit(nopython=True, cache=True)
def find_tops(coordinates, size, euler_angles):
    R = culculate_R(euler_angles)
    R = np.linalg.inv(R)
    tops = np.asarray([
        [0,         0,          0       ],
        [size[0],   0,          0       ],
        [size[0],   size[1],    0       ],
        [0,         size[1],    0       ],
        [0,         0,          size[2] ],
        [size[0],   0,          size[2] ],
        [size[0],   size[1],    size[2] ],
        [0,         size[1],    size[2] ]
    ], dtype=np.float)
    rotating_the_coordinates(tops, R)
    for i in range(tops.shape[0]):
        tops[i] += coordinates
    return tops


def inside_region(coordinates, tops):
    """
    tops  =  numpy array of the shape (8,3) with coordinates in the clockwise order. first the bottom plane is considered then the top one.
    coordinates = array of points with shape (N, 3).

    Returns the indices of the points array which are outside the cube3d
    """
    b4, b1, b2, b3, t4, t1, t2, t3 = tops

    dir1 = (t1 - b1)
    size1 = np.linalg.norm(dir1)
    dir1 = dir1/size1

    dir2 = (b2 - b1)
    size2 = np.linalg.norm(dir2)
    dir2 = dir2/size2

    dir3 = (b4 - b1)
    size3 = np.linalg.norm(dir3)
    dir3 = dir3/size3

    cube3d_center = (b1 + t3)/2.0

    dir_vec = coordinates - cube3d_center

    res1 = (np.absolute(np.dot(dir_vec, dir1))*2) <= size1
    res2 = (np.absolute(np.dot(dir_vec, dir2))*2) <= size2
    res3 = (np.absolute(np.dot(dir_vec, dir3))*2) <= size3

    indices = np.nonzero((res1*res2*res3))[0]

    return indices


@jit(nopython=True, cache=True)
def generate_directions(n):
    directions = np.empty((n, 3))
    for i in range(n):
        a1 = random()
        a2 = random()
        cos_alpha = 1 - 2*a1
        sq = sqrt(1 - cos_alpha**2)
        cos_beta = sq*cos(2*pi*a2)
        cos_gamma = sq*sin(2*pi*a2)
        directions[i] = (cos_alpha, cos_beta, cos_gamma)
    return directions


if __name__ == "__main__":
    import pyqtgraph as pg
    import pyqtgraph.opengl as gl
    from pyqtgraph.Qt import QtCore, QtGui

    directions = generate_directions(10**6)
    print('!!!')
    
    pg.mkQApp() 
    spaceView = gl.GLViewWidget()
    spaceView.setGeometry(0, 110, 1920, 1080)
    spaceView.setCameraPosition(distance=5000)
    volume = gl.GLScatterPlotItem()
    spaceView.addItem(volume)
    volume.setData(pos=directions*1000, color=(0, 0, 255, 255), size=1)

    spaceView.show()
    QtGui.QApplication.exec()

