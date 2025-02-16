import numpy as np
from numpy import load, cos, sin, log, sqrt
from particles import Photons
from hepunits import*


class Source:
    """
    Класс источника частиц

    [coordinates = (x, y, z)] = cm

    [activity] = Bq

    [distribution] = float[:,:,:]

    [voxel_size] = cm

    [energy] = eV

    [half_life] = sec
    """

    def __init__(self, coordinates, activity, distribution, voxel_size=0.4, radiation_type='Gamma', energy=140.*10**3, half_life=6*60*60, rotation_angles=None, rotation_center=None):
        self.coordinates = np.asarray(coordinates)
        self.initial_activity = np.asarray(activity)
        self.distribution = np.asarray(distribution)
        self.distribution /= np.sum(self.distribution)
        self.voxel_size = voxel_size
        self.size = np.asarray(self.distribution.shape)*self.voxel_size
        self.radiation_type = radiation_type
        self.energy = energy
        self.half_life = half_life
        self.timer = 0.
        self._rotated = False
        self._generate_emission_table()
        self.rotate(rotation_angles, rotation_center)
        self.rng = np.random.default_rng()

    def rotate(self, rotation_angles, rotation_center=None):
        if rotation_angles is not None:
            self._rotated = True
        else:
            rotation_angles = [0., 0., 0.]
        self.rotation_angles = np.asarray(rotation_angles)
        if rotation_center is None:
            rotation_center = np.asarray(self.size/2)
        self.rotation_center = rotation_center
        alpha, beta, gamma = -self.rotation_angles
        R = np.asarray([
            [cos(alpha)*cos(beta),  cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),    cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma) ],
            [sin(alpha)*cos(beta),  sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),    sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma) ],
            [-sin(beta),            cos(beta)*sin(gamma),                                       cos(beta)*cos(gamma)                                    ]
        ])
        self.R = R.T
        self._generate_emission_table()

    def _generate_emission_table(self):
        xs, ys, zs = np.meshgrid(
            np.linspace(0, self.size[0], self.distribution.shape[0], endpoint=False),
            np.linspace(0, self.size[1], self.distribution.shape[1], endpoint=False),
            np.linspace(0, self.size[2], self.distribution.shape[2], endpoint=False),
            indexing = 'ij'
        )
        coordinates = np.stack((xs, ys, zs), axis=3).reshape(-1, 3)
        probability = self.distribution.ravel()
        indices = probability.nonzero()[0]
        self.emission_table = [coordinates[indices], probability[indices]]

    @property
    def activity(self):
        return self.initial_activity*2**(-self.timer/self.half_life)
    
    @property
    def nuclei_number(self):
        return self.activity*self.half_life/log(2)

    def set_state(self, timer, rng_state=None):
        if timer is not None:
            self.timer = timer
        if rng_state is None:
            return
        self.rng.bit_generator.state['state'] = rng_state

    def generate_coordinates(self, n):
        coordinates = self.emission_table[0]
        probability = self.emission_table[1]
        coordinates = self.rng.choice(coordinates, n, p=probability)
        coordinates += self.rng.uniform(0., self.voxel_size, coordinates.shape)
        if self._rotated:
            coordinates -= self.rotation_center
            np.matmul(coordinates, self.R, out=coordinates)
            coordinates += self.rotation_center
        coordinates += self.coordinates
        return coordinates

    def generate_emission_time(self, n):
        dt = log((self.nuclei_number + n)/self.nuclei_number)*self.half_life/log(2)
        a = 2**(-self.timer/self.half_life)
        b = 2**(-(self.timer + dt)/self.half_life)
        alpha = self.rng.uniform(b, a, n)
        emission_time = -log(alpha)*self.half_life/log(2)
        return emission_time, dt

    def generate_directions(self, n):
        a1 = self.rng.random(n)
        a2 = self.rng.random(n)
        cos_alpha = 1 - 2*a1
        sq = sqrt(1 - cos_alpha**2)
        cos_beta = sq*cos(2*pi*a2)
        cos_gamma = sq*sin(2*pi*a2)
        directions = np.column_stack((cos_alpha, cos_beta, cos_gamma))
        return directions

    def generate_particles(self, n):
        energies = np.full(n, self.energy)
        directions = self.generate_directions(n)
        coordinates = self.generate_coordinates(n)
        emission_time, dt = self.generate_emission_time(n)
        self.timer += dt
        particles = Photons(energies, directions, coordinates, emission_time)
        return particles


class PointSource(Source):
    """
    Источник 99mТс-MIBI

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    
    [energy] = eV
    """

    def __init__(self, coordinates, activity, energy):
        distribution = [[[1.]]]
        voxel_size = 1.*mm
        radiation_type = 'Gamma'
        half_life = 6.*hour
        rotation_angles = None
        rotation_center = None
        super().__init__(coordinates, activity, distribution, voxel_size, radiation_type, energy, half_life, rotation_angles, rotation_center)

    def generate_coordinates(self, n):
        return np.full((n, 3), self.coordinates)

class Тс99m_MIBI(Source):
    """
    Источник 99mТс-MIBI

    [coordinates = (x, y, z)] = cm

    [activity] = Bq

    [distribution] = float[:,:,:]

    [voxel_size] = cm
    """

    def __init__(self, coordinates, activity, distribution, voxel_size, rotation_angles=None, rotation_center=None):
        radiation_type = 'Gamma'
        energy = 140.5*keV
        half_life = 6.*hour
        super().__init__(coordinates, activity, distribution, voxel_size, radiation_type, energy, half_life, rotation_angles, rotation_center)


class SourcePhantom(Тс99m_MIBI):
    """
    Источник 99mТс-MIBI

    [coordinates = (x, y, z)] = cm

    [activity] = Bq

    [phantom_name] = string

    [voxel_size] = cm
    """

    def __init__(self, coordinates, activity, phantom_name, voxel_size, rotation_angles=None, rotation_center=None):
        distribution = load(f'Phantoms/{phantom_name}.npy')
        super().__init__(coordinates, activity, distribution, voxel_size, rotation_angles=rotation_angles, rotation_center=rotation_center)


class efg3(SourcePhantom):
    """
    Источник efg3

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    """

    def __init__(self, coordinates, activity, rotation_angles=None, rotation_center=None):
        phantom_name = 'efg3'
        voxel_size = 4.*mm
        super().__init__(coordinates, activity, phantom_name, voxel_size, rotation_angles, rotation_center)
        

class efg3cut(SourcePhantom):
    """
    Источник efg3cut

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    """

    def __init__(self, coordinates, activity, rotation_angles=None, rotation_center=None):
        phantom_name = 'efg3cut'
        voxel_size = 4.*mm
        super().__init__(coordinates, activity, phantom_name, voxel_size, rotation_angles, rotation_center)


class efg3cutDefect(SourcePhantom):
    """
    Источник efg3cutDefect

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    """

    def __init__(self, coordinates, activity, rotation_angles=None, rotation_center=None):
        phantom_name = 'efg3cutDefect'
        voxel_size = 4.*mm
        super().__init__(coordinates, activity, phantom_name, voxel_size, rotation_angles, rotation_center)

