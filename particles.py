import numpy as np
from numpy import cos, sin, abs


class Particles:
    """ 
    Класс частиц

    [enegries] = eV

    [direction] = cm

    [coordinates] = cm
    """

    processes = []

    def __init__(self, energy, direction, coordinates):
        self._energy = np.asarray(energy)
        self._direction = np.asarray(direction)
        self._coordinates = np.asarray(coordinates)
        self._distance_traveled = np.zeros_like(self._energy)

    def move(self, distance):
        """ Переместить частицы """
        self._distance_traveled += distance
        self._coordinates += self._direction*distance.reshape((-1, 1))

    def rotate(self, theta, phi, indices):
        """
        Повернуть направления частиц

        [theta] = radian
        [phi] = radian
        """
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        delta1 = sin_theta*cos(phi)
        delta2 = sin_theta*sin(phi)
        delta = np.ones_like(cos_theta) - 2*(self._direction[indices, 2] < 0)
        b = self._direction[indices, 0]*delta1 + self._direction[indices, 1]*delta2
        tmp = cos_theta - b/(1 + abs(self._direction[indices, 2]))
        cos_alpha = self._direction[indices, 0]*tmp + delta1
        cos_beta = self._direction[indices, 1]*tmp + delta2
        cos_gamma = self._direction[indices, 2]*cos_theta - delta*b
        self._direction[indices] = np.column_stack((cos_alpha, cos_beta, cos_gamma))

    def change_energy(self, energy_change, indices):
        """ Измененить энергии частиц """
        self._energy[indices] -= energy_change

    def add(self, particles):
        """ Добавить частицы """
        self._energy= np.concatenate([self._energy, particles._energy])
        self._direction = np.concatenate([self._direction, particles.direction])
        self._coordinates = np.concatenate([self._coordinates, particles.coordinates])
        self._distance_traveled = np.concatenate([self._distance_traveled, particles._distance_traveled])

    def delete(self, indices):
        """ Удалить частицы """
        self._energy = np.delete(self._energy, indices)
        self._direction = np.delete(self._direction, indices, 0)
        self._coordinates = np.delete(self._coordinates, indices, 0)
        self._distance_traveled = np.delete(self._distance_traveled, indices)

    def replace(self, particles, indices):
        self._energy[indices] = particles._energy
        self._direction[indices] = particles._direction
        self._coordinates[indices] = particles.coordinates
        self._distance_traveled[indices] = particles._distance_traveled

    @property
    def energy(self):
        """ Энергии частиц """
        return self._energy

    @property
    def direction(self):
        """ Направления частиц """
        return self._direction

    @property
    def coordinates(self):
        """ Координаты частиц """
        return self._coordinates

    @property
    def distance_traveled(self):
        """ Расстояние пройденное частицами """
        return self._distance_traveled

    @property
    def count(self):
        """ Число частиц """
        return self._energy.size


class Photons(Particles):
    """ 
    Класс фотонов

    [enegries] = eV

    [direction] = cm
    
    [coordinates] = cm

    [emission_time] = sec
    """

    processes = ['PhotoelectricEffect', 'ComptonScattering']

    def __init__(self, energy, direction, coordinates, emission_time):
        super().__init__(energy, direction, coordinates)
        self._emission_time = np.asarray(emission_time)
        self._emission_coordinates = self.coordinates.copy()

    def add(self, particles):
        """ Добавить частицы """
        super().add(particles)
        self._emission_time = np.concatenate([self._emission_time, particles._emission_time])
        self._emission_coordinates = np.concatenate([self._emission_coordinates,  particles._emission_coordinates])

    def delete(self, indices):
        """ Удалить частицы """
        super().delete(indices)
        self._emission_time = np.delete(self._emission_time, indices)
        self._emission_coordinates = np.delete(self._emission_coordinates, indices, 0)

    def replace(self, particles, indices):
        """ Заменить частицы """
        super().replace(particles, indices)
        self._emission_time[indices] = particles._emission_time
        self._emission_coordinates[indices] = particles._emission_coordinates

    @property
    def emission_time(self):
        """ Времена эмиссии частиц """
        return self._emission_time

    @property
    def emission_coordinates(self):
        """ Координаты эмисии """
        return self._emission_coordinates

