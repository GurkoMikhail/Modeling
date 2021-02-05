from numpy import asarray, zeros_like, column_stack, ones_like
from numpy import cos, sin, delete, copy


class Particles:
    """ 
    Класс частиц

    [enegries] = eV

    [direction] = cm

    [coordinates] = cm
    """

    processes = []

    def __init__(self, energy, direction, coordinates):
        self._energy = asarray(energy)
        self._direction = asarray(direction)
        self._coordinates = asarray(coordinates)
        self._distance_traveled = zeros_like(self._energy)

    def move(self, distance, indices=None):
        """ Переместить частицы """
        if indices is None:
            self._coordinates += self._direction*column_stack([distance]*3)
            self._distance_traveled += distance
        else:
            self._coordinates[indices] += self._direction[indices]*column_stack([distance]*3)
            self._distance_traveled[indices] += distance

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
        delta = ones_like(cos_theta) - 2*(self._direction[indices, 2] < 0)
        b = self._direction[indices, 0]*delta1 + self._direction[indices, 1]*delta2
        tmp = cos_theta - b/(1 + abs(self._direction[indices, 2]))
        cos_alpha = self._direction[indices, 0]*tmp + delta1
        cos_beta = self._direction[indices, 1]*tmp + delta2
        cos_gamma = self._direction[indices, 2]*cos_theta - delta*b
        self._direction[indices] = column_stack((cos_alpha, cos_beta, cos_gamma))

    def change_energy(self, energy_change, indices):
        """ Измененить энергии частиц """
        self._energy[indices] -= energy_change

    def delete(self, indices):
        """ Удалить частицы из рассмотрения """
        self._energy = delete(self._energy, indices)
        self._direction = delete(self._direction, indices, 0)
        self._coordinates = delete(self._coordinates, indices, 0)
        self._distance_traveled = delete(self._distance_traveled, indices)

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
        self._emission_time = asarray(emission_time)
        self._emission_coordinates = copy(self.coordinates)

    def delete(self, indices):
        """ Удалить частицы из рассмотрения """
        super().delete(indices)
        self._emission_time = delete(self._emission_time, indices)
        self._emission_coordinates = delete(self._emission_coordinates, indices, 0)

    @property
    def emission_time(self):
        """ Времена эмиссии частиц """
        return self._emission_time

    @property
    def emission_coordinates(self):
        """ Координаты эмисии """
        return self._emission_coordinates

