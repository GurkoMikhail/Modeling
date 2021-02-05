from g4compton import generation_theta, culculate_energy_change
from numpy import array, random, nonzero, pi, sum
from materials import get_lac


class Interaction:
    """ Класс взаимодействия """

    def __init__(self, particles, space):
        self.particles = particles
        self.space = space
        self.processes = []
        for process in self.particles.processes:
            self.processes.append(processes[process](particles))
        self.data = []
        self.max_lac = get_lac(array(5), array([self.particles.energy.min(), ]), self.processes).max()
        self.rng_choose = random.default_rng()
        self.rng_free_path = random.default_rng()

    def get_lac(self):
        material = self.space.get_material(self.particles.coordinates)
        lac = get_lac(material, self.particles.energy, self.processes)
        return lac

    def get_free_path(self):
        free_path = self.rng_free_path.exponential(1/self.max_lac, self.particles.count)
        return free_path

    def choose(self, interaction_probability):
        indices = []
        rnd = self.rng_choose.random(interaction_probability[0].size)
        p0 = 0
        for p in interaction_probability:
            p1 = p0 + p
            in_delta = (p0 <= rnd)
            in_delta *= (rnd <= p1)
            ind = nonzero(in_delta)[0]
            indices.append(ind)
            p0 = p1
        return indices
        
    def apply(self):
        if self.particles.count:
            lac = self.get_lac()
            self.max_lac = sum(lac, axis=0).max()
            interaction_probability = lac/self.max_lac
            interacted = self.choose(interaction_probability)
            for i, process in enumerate(self.processes):
                indices = interacted[i]
                self.data.append(process.apply(indices))
        

class Process:
    """ Класс процесса """

    def __init__(self, particles):
        """ Конструктор процесса """
        
        self.particles = particles

    def apply(self, interacted):
        """ Применить процесс """
        name = self.__class__.__name__
        coordinates = self.particles.coordinates[interacted]
        emission_time = self.particles.emission_time[interacted]
        distance_traveled = self.particles.distance_traveled[interacted]
        emission_coordinates = self.particles.emission_coordinates[interacted]
        data = {
            'Process': name,
            'Coordinates': coordinates,
            'Emission time': emission_time,
            'Distance traveled': distance_traveled,
            'Emission coordinates': emission_coordinates
        }
        return data


class PhotoelectricEffect(Process):
    """ Класс фотоэффекта """

    def apply(self, interacted):
        """ Применить фотоэффект """
        energy_change = self.particles.energy[interacted]
        self.particles.change_energy(energy_change, interacted)
        data = super().apply(interacted)
        data.update({'Energy transfer': energy_change})
        return data
        

class ComptonScattering(Process):
    """ Класс эффекта Комптона """

    def __init__(self, particles, **kwds):
        super().__init__(particles, **kwds)
        self.rng_phi = random.default_rng()

    def get_theta(self, interacted):
        """ Получить угл рассеяния - theta """
        theta = generation_theta(self.particles.energy[interacted])
        return theta

    def get_phi(self, interacted):
        """ Получить угл рассеяния - phi """
        phi = pi*(self.rng_phi.random(interacted.size)*2 - 1)
        return phi

    def culculate_energy_change(self, theta, interacted):
        """ Вычислить изменения энергий """
        energy_change = culculate_energy_change(self.particles.energy[interacted], theta)
        return energy_change

    def apply(self, interacted):
        """ Применить эффект Комптона """
        theta = self.get_theta(interacted)
        phi = self.get_phi(interacted)
        energy_change = self.culculate_energy_change(theta, interacted)
        self.particles.rotate(theta, phi, interacted)
        self.particles.change_energy(energy_change, interacted)
        data = super().apply(interacted)
        data.update({'Energy transfer': energy_change})
        return data


class PairProduction(Process):
    """ Класс эффекта образования электрон-позитронных пар """

    pass


processes = {
    'PhotoelectricEffect': PhotoelectricEffect,
    'ComptonScattering': ComptonScattering,
    'PairProduction': PairProduction
}

