from utilites import generate_directions
import numpy as np
from particles import Photons
from processes import Interaction
import h5py
from time import time


class Modeling:
    """ 
    Основной класс моделирования
    """

    def __init__(self, space, source, **kwds):
        self.space = space
        self.source = source
        self.solid_angle = ((0, -1, 0), 10*np.pi/180)
        self.file_name = 'efg3 front projection'
        self.args = [
            'spacing',
            'solid_angle',
            'file_name'
            ]

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def start(self, total_time, time_step=1):
        for t in np.arange(self.source.timer, total_time, time_step):
            t = round(t, 4)
            dt = round(t + time_step, 4)
            flow_name = f'{(t, dt)}'
            flow = self.source.generate_particles_flow(self.space, time_step)
            flow.off_the_solid_angle(*self.solid_angle)
            print(f'Start flow for t = {t, dt}')
            start = time()
            flow.run()
            print(f'Finish flow for t = {t, dt}')
            print(f'Time left {time() - start}')
            data = self.data_structuring(flow)
            self.save_data(data, flow_name)
            self.save_data(self.source.particles_emitted, 'Source distribution')

    def data_structuring(self, flow, subject=None):
        data = {
            'Coordinates': [],
            'Energy transfer': [],
            'Emission time': [],
            'Emission coordinates': []
        }
        if subject is None:
            for dat in flow.interaction.data:
                data['Coordinates'].append(dat['Coordinates'])
                data['Energy transfer'].append(dat['Energy transfer'])
                data['Emission time'].append(dat['Emission time'])
                data['Emission coordinates'].append(dat['Emission coordinates'])
        else:
            for dat in flow.interaction.data:
                coordinates = dat['Coordinates']
                energy_transfer = dat['Energy transfer']
                emission_time = dat['Emission time']
                emission_coordinates = dat['Emission coordinates']
                inside_subject = subject.inside(coordinates)
                if inside_subject.size > 0:
                    data['Coordinates'].append(coordinates[inside_subject])
                    data['Energy transfer'].append(energy_transfer[inside_subject])
                    data['Emission time'].append(emission_time[inside_subject])
                    data['Emission coordinates'].append(emission_coordinates[inside_subject])
        data['Coordinates'] = np.concatenate(data['Coordinates'])
        data['Energy transfer'] = np.concatenate(data['Energy transfer'])
        data['Emission time'] = np.concatenate(data['Emission time'])
        data['Emission coordinates'] = np.concatenate(data['Emission coordinates'])
        return data

    def save_data(self, data, name):
        # self.lock.acquire()
        file = h5py.File(f'Output data/{self.file_name}', 'a')
        if type(data) is dict:
            group = file.create_group(f'Flows/{name}')
            for key in data.keys():
                dset = group.create_dataset(str(key), data=data[key])
        else:
            if str(name) in file:
                file[str(name)][:] = data
            else:
                dset = file.create_dataset(str(name), data=data)
        file.close()


class ParticleFlow:
    """ Класс потока частиц """

    def __init__(self, particles, space):
        self.particles = particles
        self.space = space
        self.interaction = Interaction(particles, space)
        self.left_the_space = 0
        self.step = 1
        self.min_energy = 0

    def low_energy(self):
        indices = np.nonzero(self.particles.energy <= self.min_energy)[0]
        self.particles.delete(indices)
        return indices

    def off_the_space(self):
        coordinates = self.particles.coordinates
        size = self.space.size
        x = np.nonzero((coordinates[:, 0] > size[0]) + (coordinates[:, 0] < 0))[0]
        y = np.nonzero((coordinates[:, 1] > size[1]) + (coordinates[:, 1] < 0))[0]
        z = np.nonzero((coordinates[:, 2] > size[2]) + (coordinates[:, 2] < 0))[0]
        indices = np.hstack((x, y, z))
        self.left_the_space += indices.size
        self.particles.delete(indices)
        return indices

    def off_the_solid_angle(self, vector, angle):
        cos_alpha = vector[0]*self.particles.direction[:, 0]
        cos_alpha += vector[1]*self.particles.direction[:, 1]
        cos_alpha += vector[2]*self.particles.direction[:, 2]
        # cos_alpha = x + y + z
        indices = np.nonzero(cos_alpha <= np.cos(angle))[0]
        self.particles.delete(indices)
        return indices

    def next_step(self):
        free_path = self.interaction.get_free_path()
        self.particles.move(free_path)
        self.off_the_space()
        self.interaction.apply()
        self.low_energy()
        self.step += 1
    
    def run(self):
        """ Реализация работы процесса """
        while self.particles.count:
                self.next_step()


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

    def __init__(self, coordinates, activity, distribution, voxel_size=0.4, radiation_type='Gamma', energy=140.*10**3, half_life=6*60*60, **kwds):
        self.coordinates = np.asarray(coordinates)
        self.initial_activity = np.asarray(activity)
        self.distribution = np.asarray(distribution)
        self.distribution /= np.sum(self.distribution)
        self.particles_emitted = np.zeros_like(self.distribution)
        self.voxel_size = voxel_size
        self.radiation_type = radiation_type
        self.energy = energy
        self.half_life = half_life
        self.timer = 0
        self.coordinates_table = self.generate_coordinates_table()
        self.args = []
        self.rng_dist = np.random.default_rng()
        self.rng_time = np.random.default_rng()
        self.rng_dir = np.random.default_rng()

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def generate_coordinates_table(self):
        coordinates_table = []
        size = np.asarray(self.distribution.shape, dtype=np.float)
        size *= self.voxel_size
        size += self.coordinates
        for x in np.arange(self.coordinates[0], size[0], self.voxel_size):
            for y in np.arange(self.coordinates[1], size[1], self.voxel_size):
                for z in np.arange(self.coordinates[2], size[2], self.voxel_size):
                    coordinates_table.append([x, y, z])
        coordinates_table = np.asarray(coordinates_table)
        return coordinates_table

    @property
    def activity(self):
        return self.initial_activity*2**(-self.timer/self.half_life)
    
    @property
    def nuclei_number(self):
        return self.activity*self.half_life/np.log(2)

    def save_particles_emmited(self, coordinates):
        coordinates = coordinates - self.coordinates
        coordinates /= self.voxel_size
        coordinates = coordinates.astype(np.uint64)
        self.particles_emitted[(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])] += 1

    def generate_coordinates(self, n):
        p = self.distribution.ravel()
        coordinates = self.rng_dist.choice(self.coordinates_table, n, p=p)
        return coordinates

    def generate_emission_time(self, n):
        emission_time = self.rng_time.random(n) + self.timer
        return emission_time

    def generate_directions(self, n):
        a1 = self.rng_dir.random(n)
        a2 = self.rng_dir.random(n)
        cos_alpha = 1 - 2*a1
        sq = np.sqrt(1 - cos_alpha**2)
        cos_beta = sq*np.cos(2*np.pi*a2)
        cos_gamma = sq*np.sin(2*np.pi*a2)
        directions = np.column_stack((cos_alpha, cos_beta, cos_gamma))
        return directions

    def generate_particles(self, n):
        energies = np.full(n, self.energy)
        directions = generate_directions(n)
        coordinates = self.generate_coordinates(n)
        emission_time = self.generate_emission_time(n)
        particles = Photons(energies, directions, coordinates, emission_time)
        return particles

    def generate_particles_flow(self, space, time_step):
        n = int(self.nuclei_number*(1 - 2**(-time_step/self.half_life)))
        particles = self.generate_particles(n)
        self.save_particles_emmited(particles.coordinates)
        particles_flow = ParticleFlow(particles, space)
        self.timer += time_step
        return particles_flow
