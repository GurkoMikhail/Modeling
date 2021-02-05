from numba import jit
from numpy import array, cos, sin, sqrt, dot, empty, pi
from random import random


def culculate_R(euler_angles):
    alpha, beta, gamma = euler_angles

    Rx = array([
        [1,                 0,              0               ],
        [0,                 cos(beta),   sin(beta)    ],
        [0,                 -sin(beta),  cos(beta)    ]
    ])

    Ry = array([
        [cos(gamma),     0,              -sin(gamma)  ],
        [0,                 1,              0               ],
        [sin(gamma),     0,              cos(gamma)   ]
    ])

    Rz = array([
        [cos(alpha),     sin(alpha),  0               ],
        [-sin(alpha),    cos(alpha),  0               ],
        [0,                 0,              1               ]
    ])

    Rzx = dot(Rz, Rx)
    R = dot(Rzx, Ry)
    return R


@jit(nopython=True, cache=True)
def rotating_the_coordinates(coordinates, R):
    for i in range(coordinates.shape[0]):
        coordinates[i] = dot(R, coordinates[i])


@jit(nopython=True, cache=True)
def generate_directions(n):
    directions = empty((n, 3))
    for i in range(n):
        a1 = random()
        a2 = random()
        cos_alpha = 1 - 2*a1
        sq = sqrt(1 - cos_alpha**2)
        cos_beta = sq*cos(2*pi*a2)
        cos_gamma = sq*sin(2*pi*a2)
        directions[i] = (cos_alpha, cos_beta, cos_gamma)
    return directions

