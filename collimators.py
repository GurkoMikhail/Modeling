from numba import jit
import numpy as np
from numpy import ubyte, sqrt, cos, pi


@jit(nopython=True, cache=True)
def get_collimated(coordinates, hole_centers, hole_diameter):
    collimated = np.zeros(coordinates.shape[0], dtype=ubyte)
    find_collimated(collimated, coordinates, hole_centers, hole_diameter)
    return np.nonzero(collimated)[0]


@jit(nopython=True, cache=True)
def find_collimated(collimated, coordinates, hole_centers, hole_diameter):
    a = sqrt(3)/4
    a_2 = a/2
    for i in range(coordinates.shape[0]):
        x0 = coordinates[i, 0]
        y0 = coordinates[i, 1]
        for hole_center in hole_centers:
            x = hole_center[0]
            y = hole_center[1]
            dx = abs(x - x0)/hole_diameter
            dy = abs(y - y0)/hole_diameter
            if dy <= a and a*dx + dy/4 <= a_2:
                collimated[i] = True
                break


@jit(nopython=True, cache=True)
def get_centers(size, hole_diameter, septa):
    hole_size_X = hole_diameter/2
    hole_size_Y = (hole_size_X*sqrt(3)/2)
    intercenter_distance_X = (hole_size_X*(1 + cos(pi/3)) + septa)*2
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
    centers = np.array(centers)
    dx = (size[0] - centers[-1, 0] - (hole_size_X + septa))
    dy = (size[1] - centers[-1, 1] - (hole_size_Y + septa))
    centers[:, 1] += dy/2
    if dx > hole_size_X + septa:
        dx -= intercenter_distance_X/2
    centers[:, 0] += dx/2
    return centers
