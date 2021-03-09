import g4compton
import numpy as np
from numpy import pi
import materials


class Interaction:
    """ Класс взаимодействия """

    def __init__(self, particles, space):
        self.particles = particles
        self.space = space
        self.processes = []
        for process in self.particles.processes:
            self.processes.append(processes[process](particles))
        self.data = []
        self.space_materials_list = self.space.materials_list
        self.rng_choose = np.random.RandomState()
        self.rng_free_path = np.random.RandomState()

    def get_lac(self, indices):
        coordinates = self.particles.coordinates[indices]
        energy = self.particles.energy[indices]
        material = self.space.get_material(coordinates)
        lac = materials.get_lac(material, energy, self.processes)
        return lac

    def get_free_path(self):
        self.max_lac = materials.get_max_lac(self.space_materials_list, self.particles.energy, self.processes)
        free_path = self.rng_free_path.exponential(1/self.max_lac, self.particles.count)
        return free_path

    def choose(self, selectable, probability):
        indices = []
        rnd = self.rng_choose.rand(probability[0].size)
        p0 = 0
        for p in probability:
            p1 = p0 + p
            in_delta = (p0 <= rnd)
            in_delta *= (rnd <= p1)
            ind = np.nonzero(in_delta)[0]
            indices.append(selectable[ind])
            p0 = p1
        return indices
        
    def apply(self, indices):
        if indices.size:
            lac = self.get_lac(indices)
            max_lac = self.max_lac[indices]
            interaction_probability = lac/max_lac
            interacted = self.choose(indices, interaction_probability)
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
        

class CoherentScattering(Process):
    """ Класс когерентного рассеяния """
    
    def __init__(self, particles, **kwds):
        super().__init__(particles, **kwds)
        self.rng_phi = np.random.RandomState()


class ComptonScattering(Process):
    """ Класс эффекта Комптона """

    def __init__(self, particles, **kwds):
        super().__init__(particles, **kwds)
        self.rng_phi = np.random.RandomState()

    def get_theta(self, interacted):
        """ Получить угл рассеяния - theta """
        theta = g4compton.generation_theta(self.particles.energy[interacted])
        return theta

    def get_phi(self, interacted):
        """ Получить угл рассеяния - phi """
        phi = pi*(self.rng_phi.rand(interacted.size)*2 - 1)
        return phi

    def culculate_energy_change(self, theta, interacted):
        """ Вычислить изменения энергий """
        energy_change = g4compton.culculate_energy_change(self.particles.energy[interacted], theta)
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

