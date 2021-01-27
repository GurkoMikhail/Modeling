import numpy as np
import materials
from time import time

class Interaction:
    """ Класс взаимодействия """

    def __init__(self, particles, space, **kwds):
        self.particles = particles
        self.space = space
        self.processes = []
        for process in self.particles.processes:
            self.processes.append(processes[process](particles))
        self.data = []
        self.rng_choose = np.random.default_rng()
        self.rng_free_path = np.random.default_rng()
        self.args = []

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def probability(self, travel_distance):
        material = self.space.get_material(self.particles.coordinates)
        lac = materials.get_lac(material, self.particles.energy, self.processes)
        probability = 1 - np.exp(-lac*travel_distance)
        return probability

    def choose(self, interaction_probability):
        indices = []
        rnd = self.rng_choose.random(interaction_probability[0].size)
        p0 = 0
        for p in interaction_probability:
            p1 = p0 + p
            in_delta = (p0 <= rnd)
            in_delta *= (rnd <= p1)
            ind = np.nonzero(in_delta)[0]
            indices.append(ind)
            p0 = p1
        return indices

    def interacted_list(self, travel_distance):
        interaction_probability = self.probability(travel_distance)
        interacted = self.choose(interaction_probability)
        for indices in interacted:
            travel_distance[indices] *= self.rng_free_path.random(indices.size)
        return interacted

    def apply(self, interacted):
        for i, process in enumerate(self.processes):
            indices = interacted[i]
            self.data.append(process.apply(indices))
        

class Process:
    """ Класс процесса """

    def __init__(self, particles, **kwds):
        """ Конструктор процесса """
        
        self.particles = particles
        self.args = []

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

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


import photoelectriceffect 

class PhotoelectricEffect(Process):
    """ Класс фотоэффекта """

    def apply(self, interacted):
        """ Применить фотоэффект """
        # print(f'PhotoelectricEffect {interacted.size}')
        energy_change = self.particles.energy[interacted]
        self.particles.change_energy(energy_change, interacted)
        data = super().apply(interacted)
        data.update({'Energy transfer': energy_change})
        return data


import g4compton

class ComptonScattering(Process):
    """ Класс эффекта Комптона """

    def get_theta(self, interacted):
        """ Получить угл рассеяния - theta """
        theta = g4compton.generation_theta(self.particles.energy[interacted])
        return theta

    def get_phi(self, interacted):
        """ Получить угл рассеяния - phi """
        phi = np.pi*(np.random.random_sample(interacted.size)*2 - 1)
        return phi

    def culculate_energy_change(self, theta, interacted):
        """ Вычислить изменения энергий """
        energy_change = g4compton.culculate_energy_change(self.particles.energy[interacted], theta)
        return energy_change

    def apply(self, interacted):
        """ Применить эффект Комптона """
        # print(f'Compton {interacted.size}')
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

