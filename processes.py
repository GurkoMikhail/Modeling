import g4compton
import g4coherent
import numpy as np
from numpy import pi, cos


class Interaction:
    """ Класс взаимодействия """

    def __init__(self, particles, space, materials):
        self.particles = particles
        self.space = space
        self.materials = materials
        self.processes = []
        for process in self.particles.processes:
            self.processes.append(processes[process](self.particles, self.materials))
        lac_funtions = self.materials.construct_lac_funtions(self.processes)
        self.lacs = lac_funtions['Total']
        for process in self.processes:
            process.lacs = lac_funtions[process.name]
        self.rng_choose = np.random.default_rng()
        self.rng_free_path = np.random.default_rng()

    def casting(self):
        subject_index, distance = self.space.ray_casting(
            self.particles.coordinates,
            self.particles.direction
            )
        material, complex_subject = self.space.get_material_of_subject(subject_index)
        total_lac = self.get_total_lac(material, self.particles.energy)
        free_path = self.get_free_path(total_lac)
        interacted = (free_path < distance).nonzero()[0]
        distance[interacted] = free_path[interacted]
        self.particles.move(distance)
        inside_space = self.space.inside(self.particles.coordinates[interacted])
        interacted = interacted[inside_space]
        interaction_data = {}
        if interacted.size:
            coordinates = self.particles.coordinates[complex_subject]
            material[complex_subject] = self.space.get_material(coordinates)
            material = material[interacted]
            total_lac = total_lac[interacted]
            interaction_data.update(self.apply(interacted, material, total_lac))
        return interaction_data

    def get_processes_lac(self, materials, energy):
        lac_out = np.zeros((len(self.processes), energy.size))
        for material in np.unique(materials):
            indices = (materials == material).nonzero()[0]
            for i, process in enumerate(self.processes):
                lac_out[i, indices] = process.lacs[material](energy[indices])
        return lac_out

    def get_total_lac(self, materials, energy):
        total_lac = np.empty_like(energy)
        for material in np.unique(materials):
            indices = (materials == material).nonzero()[0]
            total_lac[indices] = self.lacs[material](energy[indices])
        return total_lac

    def get_free_path(self, total_lac):
        free_path = self.rng_free_path.exponential(1/total_lac, total_lac.size)
        return free_path

    def choose_interaction(self, probabilities):
        indices = []
        rnd = self.rng_choose.random(probabilities[0].size)
        p0 = 0
        for probability in probabilities:
            p1 = p0 + probability
            in_delta = (p0 <= rnd)
            in_delta *= (rnd <= p1)
            ind = in_delta.nonzero()[0]
            indices.append(ind)
            p0 = p1
        return indices
    
    def apply(self, indices, material, total_lac):
        data = {}
        energy = self.particles.energy[indices]
        lacs = self.get_processes_lac(material, energy)
        interaction_probabilities = lacs/total_lac
        interacted = self.choose_interaction(interaction_probabilities)
        for i, process in enumerate(self.processes):
            ind = interacted[i]
            process_data = process.apply(indices[ind], material[ind])
            data.update({process.name: process_data})
        return data


class Process:
    """ Класс процесса """

    def __init__(self, particles, materials):
        """ Конструктор процесса """
        self.particles = particles
        self.materials = materials

    @property
    def name(self):
        return self.__class__.__name__

    def apply(self, interacted):
        """ Применить процесс """
        coordinates = self.particles.coordinates[interacted]
        emission_time = self.particles.emission_time[interacted]
        distance_traveled = self.particles.distance_traveled[interacted]
        emission_coordinates = self.particles.emission_coordinates[interacted]
        data = {
            'Coordinates': coordinates,
            'Emission time': emission_time,
            'Distance traveled': distance_traveled,
            'Emission coordinates': emission_coordinates
        }
        return data


class PhotoelectricEffect(Process):
    """ Класс фотоэффекта """

    def apply(self, interacted, materials):
        """ Применить фотоэффект """
        energy_change = self.particles.energy[interacted]
        self.particles.change_energy(energy_change, interacted)
        data = super().apply(interacted)
        data.update({'Energy transfer': energy_change})
        return data
        

class CoherentScattering(Process):
    """ Класс когерентного рассеяния """
    
    def __init__(self, particles, materials):
        super().__init__(particles, materials)
        self.rng_phi = np.random.default_rng()

    def get_theta(self, interacted, materials):
        energy = self.particles.energy[interacted]
        Z = self.materials.select_atom(materials)
        theta = g4coherent.generation_theta(energy, Z)
        return theta

    def get_phi(self, interacted):
        """ Получить угол рассеяния - phi """
        phi = pi*(self.rng_phi.random(interacted.size)*2 - 1)
        return phi

    def apply(self, interacted, materials):
        """ Применить эффект Комптона """
        theta = self.get_theta(interacted, materials)
        phi = self.get_phi(interacted)
        self.particles.rotate(theta, phi, interacted)
        data = super().apply(interacted)
        data.update({'Energy transfer': np.zeros_like(phi)})
        return data


class ComptonScattering(CoherentScattering):
    """ Класс эффекта Комптона """

    def get_theta(self, interacted):
        """ Получить угол рассеяния - theta """
        theta = g4compton.generation_theta(self.particles.energy[interacted])
        return theta

    def culculate_energy_change(self, theta, interacted):
        """ Вычислить изменения энергий """
        energy = self.particles.energy[interacted]
        k = energy/mec2
        k1_cos = k*(1 - cos(theta))
        energy_change = energy*k1_cos/(1 + k1_cos)
        return energy_change

    def apply(self, interacted, materials):
        """ Применить эффект Комптона """
        theta = self.get_theta(interacted)
        phi = self.get_phi(interacted)
        energy_change = self.culculate_energy_change(theta, interacted)
        self.particles.rotate(theta, phi, interacted)
        self.particles.change_energy(energy_change, interacted)
        data = Process.apply(self, interacted)
        data.update({'Energy transfer': energy_change})
        return data


class PairProduction(Process):
    """ Класс эффекта образования электрон-позитронных пар """

    pass


mec2 = 510998.9461          #eV

processes = {
    'PhotoelectricEffect': PhotoelectricEffect,
    'ComptonScattering': ComptonScattering,
    'CoherentScattering': CoherentScattering,
    'PairProduction': PairProduction
}

