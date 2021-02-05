from numba import vectorize, float64
from numpy import log, exp, sqrt, arccos, cos
from random import random


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
    theta = arccos(costheta)
    return theta

@vectorize([float64(float64, float64)], nopython=True, cache=True)
def culculate_energy_change(energy, theta):
    mec2 = 510998.9461          #eV
    k = energy/mec2
    k1_cos = k*(1 - cos(theta))
    energy_change = energy*k1_cos/(1+k1_cos)
    return energy_change

