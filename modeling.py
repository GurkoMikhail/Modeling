from h5py import File
import numpy as np
from numpy import pi, sqrt, cos, sin, log
import utilites
from particles import Photons
from processes import Interaction
from materials import materials_reverse_list
from time import time


class Modeling:
    """ 
    Основной класс моделирования
    """

    def __init__(self, space, source, **kwds):
        self.space = space
        self.source = source
        self.solid_angle = ((0, -1, 0), 10*pi/180)
        self.time_step = 1
        self.file_name = f'{self}'
        self.subject = None
        self.distibution_voxel_size = 0.4
        self.args = [
            'solid_angle',
            'time_step',
            'file_name',
            'subject',
            'distibution_voxel_size'
            ]

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def start(self, total_time):
        self.save_modeling_parameters()
        for t in np.arange(0, total_time, self.time_step):
            t = round(t, 5)
            dt = round(t + self.time_step, 5)
            flow_name = f'{(t, dt)}'
            if self.check_flow_in_file(flow_name):
                self.source.timer += self.time_step
            else:
                flow = self.source.generate_particles_flow(self.space, self.time_step, self.solid_angle)
                flow.run()
                self.save_data(flow)

    def check_flow_in_file(self, flow_name):
        file = File(f'Output data/{self.file_name}', 'a')
        if 'Flows' in file:
            flows = file['Flows']
            return flow_name in flows
        return False

    def save_data(self, flow):
        self.save_flow_data(flow)
        if self.subject is not None:
            self.save_dose_distribution(flow)

    def save_flow_data(self, flow):
        data = {
            'Coordinates': [],
            'Energy transfer': [],
            'Emission time': [],
            'Emission coordinates': []
        }
        if self.subject is None:
            for dat in flow.interaction.data:
                data['Coordinates'].append(dat['Coordinates'])
                data['Energy transfer'].append(dat['Energy transfer'])
                data['Emission time'].append(dat['Emission time'])
                data['Emission coordinates'].append(dat['Emission coordinates'])
        else:
            for dat in flow.interaction.data:
                coordinates = np.copy(dat['Coordinates'])
                energy_transfer = dat['Energy transfer']
                emission_time = dat['Emission time']
                emission_coordinates = dat['Emission coordinates']
                inside_subject = self.subject.inside(coordinates)
                if inside_subject.size > 0:
                    data['Coordinates'].append(coordinates[inside_subject])
                    data['Energy transfer'].append(energy_transfer[inside_subject])
                    data['Emission time'].append(emission_time[inside_subject])
                    data['Emission coordinates'].append(emission_coordinates[inside_subject])
        data['Coordinates'] = np.concatenate(data['Coordinates'])
        data['Energy transfer'] = np.concatenate(data['Energy transfer'])
        data['Emission time'] = np.concatenate(data['Emission time'])
        data['Emission coordinates'] = np.concatenate(data['Emission coordinates'])
        
        file = File(f'Output data/{self.file_name}', 'a')
        try:
            group = file.create_group(f'Flows/{flow.name}')
        except Exception:
            group = file['Flows']
        for key in data.keys():
            group.create_dataset(str(key), data=data[key])

    def save_dose_distribution(self, flow):
        file = File(f'Output data/{self.file_name}', 'a')
        try:
            group = file.create_group('Dose distribution')
        except Exception:
            group = file['Dose distribution']
            volume = group['Volume']
        else:
            volume = group.create_dataset('Volume', data=np.zeros((self.space.size/self.distibution_voxel_size).astype(np.uint), dtype=np.float64))
            group.create_dataset('Voxel size', data=self.distibution_voxel_size)
        coordinates = []
        energy_transfer = []
        for dat in flow.interaction.data:
                coordinates.append(dat['Coordinates'])
                energy_transfer.append(dat['Energy transfer'])
        coordinates = np.concatenate(coordinates)
        energy_transfer = np.concatenate(energy_transfer)
        flow_volume = np.histogramdd(
            sample=coordinates,
            bins=(self.space.size/self.distibution_voxel_size).astype(np.int),
            range=((0, self.space.size[0]), (0, self.space.size[1]), (0, self.space.size[2])),
            weights=energy_transfer
        )[0]
        volume[:] += flow_volume
        file.close()

    def save_source_distribution(self, flow):
        file = File(f'Output data/{self.file_name}', 'a')
        try:
            group = file.create_group('Source distribution')
        except Exception:
            group = file['Source distribution']
            volume = group['Volume']
        else:
            volume = group.create_dataset('Volume', data=np.zeros((self.space.size/self.distibution_voxel_size).astype(np.uint)))
            group.create_dataset('Voxel size', data=self.distibution_voxel_size)
        flow_volume = 0
        volume += flow_volume
        file.close()

    def save_modeling_parameters(self):
        file = File(f'Output data/{self.file_name}', 'a')
        try:
            group = file.create_group('Modeling parameters')
        except Exception:
            print('Параметры уже записаны')
        else:
            solidAngle = group.create_group('Solid angle')
            solidAngle.create_dataset('Vector', data=self.solid_angle[0])
            solidAngle.create_dataset('Angle', data=self.solid_angle[1])

            spaceParameters = group.create_group('Space')
            spaceParameters.create_dataset('size', data=self.space.size)
            spaceParameters.create_dataset('material', data=materials_reverse_list[self.space.material])
            for subject in self.space.subjects:
                subjectParameters = spaceParameters.create_group(subject.__class__.__name__)
                for parameter_name, value in subject.__dict__.items():
                    subjectParameters.create_dataset(parameter_name, data=value)
            
            sourceParameters = group.create_group('Source')
            for parameter_name, value in self.source.__dict__.items():
                if parameter_name not in ('particles_emitted', 'timer', 'coordinates_table', 'size', 'rng_dist', 'rng_time', 'rng_dir', 'rng_ddist'):
                    sourceParameters.create_dataset(parameter_name, data=value)

            if subject is not None:
                subject = self.subject.__class__.__name__
            group.create_dataset('Subject', data=subject)
            group.create_dataset('Interactions', data=Photons.processes)
        finally:
            file.close()


class ParticleFlow:
    """ Класс потока частиц """

    def __init__(self, particles, space, solid_angle, name):
        self.particles = particles
        self.space = space
        self.solid_angle = solid_angle
        self.name = name
        self.interaction = Interaction(particles, space)
        self.left_the_space = 0
        self.step = 1
        self.min_energy = 0

    def low_energy(self):
        indices = np.nonzero(self.particles.energy <= self.min_energy)[0]
        return indices
        
    def off_the_solid_angle(self):
        if self.solid_angle is None:
            return []
        vector, angle = self.solid_angle
        cos_alpha = vector[0]*self.particles.direction[:, 0]
        cos_alpha += vector[1]*self.particles.direction[:, 1]
        cos_alpha += vector[2]*self.particles.direction[:, 2]
        indices = np.nonzero(cos_alpha <= cos(angle))[0]
        self.particles.delete(indices)
        return indices

    @property
    def invalid_particles(self):
        indices = []
        indices.append(self.low_energy())
        indices.append(self.space.outside(self.particles.coordinates))
        indices = np.concatenate(indices)
        return indices

    @property
    def valid_particles(self):
        indices = self.space.inside(self.particles.coordinates)
        return indices

    def next_step(self):
        free_path = self.interaction.get_free_path()
        self.particles.move(free_path)
        self.interaction.apply(self.valid_particles)
        self.particles.delete(self.invalid_particles)
        self.step += 1

    def run(self):
        """ Реализация работы процесса """
        self.off_the_solid_angle()
        print(f'Start flow {self.name}')
        start = time()
        while self.particles.count:
                self.next_step()
        print(f'Finish flow flow {self.name}')
        print(f'Time left {time() - start}')


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

    def __init__(self, coordinates, activity, distribution, voxel_size=0.4, radiation_type='Gamma', energy=140.*10**3, half_life=6*60*60, euler_angles=None, rotation_center=None):
        self.coordinates = np.asarray(coordinates)
        self.initial_activity = np.asarray(activity)
        self.distribution = np.asarray(distribution)
        self.distribution /= np.sum(self.distribution)
        self.voxel_size = voxel_size
        self.size = np.asarray(self.distribution.shape)*self.voxel_size
        self.radiation_type = radiation_type
        self.energy = energy
        self.half_life = half_life
        self.timer = 0
        self.coordinates_table = self.generate_coordinates_table()
        self.rotated = False
        if euler_angles is not None:
            self.rotate(euler_angles, rotation_center)
        self.rng_dist = np.random.default_rng()
        self.rng_ddist = np.random.default_rng()
        self.rng_time = np.random.default_rng()
        self.rng_dir = np.random.default_rng()

    def rotate(self, euler_angles, rotation_center=None):
        self.rotated = True
        self.euler_angles = np.asarray(euler_angles)
        if rotation_center is None:
            rotation_center = np.asarray(self.size/2)
        self.rotation_center = rotation_center
        self.R = np.asarray(utilites.culculate_R_euler(-self.euler_angles))

    def generate_coordinates_table(self):
        coordinates_table = []
        for x in np.linspace(0, self.size[0], self.distribution.shape[0]):
            for y in np.linspace(0, self.size[1], self.distribution.shape[1]):
                for z in np.linspace(0, self.size[2], self.distribution.shape[2]):
                    coordinates_table.append([x, y, z])
        coordinates_table = np.asarray(coordinates_table)
        return coordinates_table

    @property
    def activity(self):
        return self.initial_activity*2**(-self.timer/self.half_life)
    
    @property
    def nuclei_number(self):
        return self.activity*self.half_life/log(2)

    def generate_coordinates(self, n):
        p = self.distribution.ravel()
        coordinates = self.rng_dist.choice(self.coordinates_table, n, p=p)
        dcoordinates = self.rng_ddist.uniform(0, self.voxel_size, coordinates.shape)
        coordinates += dcoordinates
        if self.rotated:
            coordinates -= self.rotation_center
            utilites.rotating_the_coordinates(coordinates, self.R)
            coordinates += self.rotation_center
        coordinates += self.coordinates
        return coordinates

    def generate_emission_time(self, n):
        emission_time = self.rng_time.random(n) + self.timer
        return emission_time

    def generate_directions(self, n):
        a1 = self.rng_dir.random(n)
        a2 = self.rng_dir.random(n)
        cos_alpha = 1 - 2*a1
        sq = sqrt(1 - cos_alpha**2)
        cos_beta = sq*cos(2*pi*a2)
        cos_gamma = sq*sin(2*pi*a2)
        directions = np.column_stack((cos_alpha, cos_beta, cos_gamma))
        return directions

    def generate_particles(self, n):
        energies = np.full(n, self.energy)
        directions = utilites.generate_directions(n)
        coordinates = self.generate_coordinates(n)
        emission_time = self.generate_emission_time(n)
        particles = Photons(energies, directions, coordinates, emission_time)
        return particles

    def generate_particles_flow(self, space, time_step, solid_angle, name=None):
        n = int(self.nuclei_number*(1 - 2**(-time_step/self.half_life)))
        particles = self.generate_particles(n)
        if name is None:
            name = f'{(round(self.timer, 5), round(self.timer + time_step, 5))}'
        particles_flow = ParticleFlow(particles, space, solid_angle, name)
        self.timer += time_step
        return particles_flow
