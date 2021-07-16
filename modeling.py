from h5py import File
import numpy as np
from numpy import pi, sqrt, cos, sin, log
import utilites
from particles import Photons
from processes import Interaction
from multiprocessing import Process, Queue
from time import time


class Modeling(Process):
    """ 
    Основной класс моделирования
    """

    def __init__(self, space, source, materials, stop_time, **kwds):
        super().__init__()
        self.space = space
        self.source = source
        self.materials = materials
        self.stop_time = stop_time
        self.start_time = 0.
        self.solid_angle = ((0, -1, 0), 20*pi/180)
        self.particles_number = 10**6
        self.flow_number = 1
        self.file_name = f'{self}'
        self.subject = None
        self.save_flows = True
        self.save_dose_data = True
        self.distibution_voxel_size = 0.4
        self.args = [
            'start_time',
            'solid_angle',
            'particles_number',
            'flow_number',
            'file_name',
            'subject',
            'save_flows',
            'save_dose_data',
            'distibution_voxel_size'
            ]

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def generate_flows(self, start_time, stop_time, flow_number=1):
        last_time = self.check_progress_in_file()
        start_time = max([start_time, last_time])
        self.source.set_timer(start_time)
        self.particles_number //= flow_number
        self.source.initial_activity //= flow_number
        queue = Queue(maxsize=flow_number)
        flows = []
        for _ in range(flow_number):
            flow = ParticleFlow(
                source=self.source,
                space=self.space,
                materials=self.materials,
                stop_time=stop_time,
                particles_number=self.particles_number,
                queue=queue,
                solid_angle=self.solid_angle
                )
            flows.append(flow)
        return flows, queue

    def run(self):
        self.source.set_timer(self.start_time)
        self.save_modeling_parameters()
        flows, queue = self.generate_flows(self.start_time, self.stop_time, self.flow_number)
        for flow in flows:
            flow.start()
        start_time = time()
        finished_flows = 0
        for step_data in iter(queue.get, None):
            if step_data == 'Finish':
                finished_flows += 1
                if finished_flows == self.flow_number:
                    print('Modeling end!')
                    break
                continue
            finish_time = time() - start_time
            self.update_step_data(step_data)
            if self.save_dose_data:
                self.update_dose_data(step_data)
            print(f'\tReal time passed: {finish_time} seconds')
            start_time = time()

    def check_progress_in_file(self):
        try:
            file = File(f'Output data/{self.file_name}', 'r')
            last_time = file['Source timer']
            last_time = float(np.array(last_time))
            file.close()
        except Exception:
            print(f'Не удалось проверить прогресс')
            last_time = 0
        finally:
            print(f'Source timer: {last_time}')
            return last_time

    def concatenate_step_data(self, step_data):
        data = {
            'Coordinates': [],
            'Energy transfer': [],
            'Emission time': [],
            'Emission coordinates': []
        }
        if self.subject is None:
            for dat in step_data['Interaction'].values():
                data['Coordinates'].append(dat['Coordinates'])
                data['Energy transfer'].append(dat['Energy transfer'])
                data['Emission time'].append(dat['Emission time'])
                data['Emission coordinates'].append(dat['Emission coordinates'])
        else:
            for dat in step_data['Interaction'].values():
                coordinates = dat['Coordinates'].copy()
                energy_transfer = dat['Energy transfer']
                emission_time = dat['Emission time']
                emission_coordinates = dat['Emission coordinates']
                inside_subject = self.subject.inside(coordinates)
                if inside_subject.size > 0:
                    data['Coordinates'].append(coordinates[inside_subject])
                    data['Energy transfer'].append(energy_transfer[inside_subject])
                    data['Emission time'].append(emission_time[inside_subject])
                    data['Emission coordinates'].append(emission_coordinates[inside_subject])
        if len(data['Coordinates']) > 0:
            data['Coordinates'] = np.concatenate(data['Coordinates'])
            data['Energy transfer'] = np.concatenate(data['Energy transfer'])
            data['Emission time'] = np.concatenate(data['Emission time'])
            data['Emission coordinates'] = np.concatenate(data['Emission coordinates'])
        return_data = {'Source timer': step_data['Source timer']}
        return_data.update({'Interactions data': data})
        return return_data

    def update_step_data(self, step_data):
        data = self.concatenate_step_data(step_data)
        source_timer = data['Source timer']
        data = data['Interactions data']
        events_number = len(data['Coordinates'])
        print(f'Source timer: {source_timer} seconds\n'
            + f'\tGenerated {events_number} events'
            )
        if events_number > 0:
            try:
                file = File(f'Output data/{self.file_name}', 'r+')
            except Exception:
                print(f'Не удалось сохранить Step №{step_data["Step"]}')
            else:
                if not 'Interactions data' in file:
                    file.create_dataset('Source timer', data=source_timer)
                    group = file.create_group('Interactions data')
                    for key in data.keys():
                        maxshape = list(data[key].shape)
                        maxshape[0] = None
                        group.create_dataset(
                            str(key),
                            data=data[key],
                            compression="gzip",
                            chunks=True,
                            maxshape=maxshape
                            )
                else:
                    timer = file['Source timer']
                    timer[...] = max([float(np.array(source_timer)), float(np.array(timer))])
                    group = file['Interactions data']
                    for key in group.keys():
                        group[key].resize(
                            (group[key].shape[0] + data[key].shape[0]),
                            axis=0
                            )
                        group[key][-data[key].shape[0]:] = data[key]

    def update_dose_data(self, step_data):
        coordinates = []
        energy_transfer = []
        for dat in step_data['Interaction'].values():
                coordinates.append(dat['Coordinates'])
                energy_transfer.append(dat['Energy transfer'])
        if len(coordinates) > 0:
            coordinates = np.concatenate(coordinates)
            energy_transfer = np.concatenate(energy_transfer)
            flow_volume = np.histogramdd(
                sample=coordinates,
                bins=(self.space.size/self.distibution_voxel_size).astype(np.int),
                range=((0, self.space.size[0]), (0, self.space.size[1]), (0, self.space.size[2])),
                weights=energy_transfer
            )[0]
            try:
                file = File(f'Output data/{self.file_name}', 'r+')
            except Exception:
                print(f'Не удалось сохранить dose Step №{step_data["Step"]}')
            else:
                if not 'Dose distribution' in file:
                    group = file.create_group('Dose distribution')
                    volume = group.create_dataset('Volume', data=np.zeros((self.space.size/self.distibution_voxel_size).astype(np.uint), dtype=np.float64))
                    group.create_dataset('Voxel size', data=self.distibution_voxel_size)
                else:
                    group = file['Dose distribution']
                    volume = group['Volume']
                volume[...] += flow_volume
                file.close()

    def save_modeling_parameters(self):
        try:
            file = File(f'Output data/{self.file_name}', 'x')
            group = file.create_group('Modeling parameters')
        except Exception:
            print(f'Не удалось записать параметры {self.file_name}')
        else:
            solidAngle = group.create_group('Solid angle')
            if self.solid_angle is not None:
                solidAngle.create_dataset('Vector', data=self.solid_angle[0])
                solidAngle.create_dataset('Angle', data=self.solid_angle[1])
            materialsGroup = group.create_group('Dict of materials indices')
            for name, index in self.materials.indices_dict.items():
                materialsGroup.create_dataset(name, data=index)
            spaceParameters = group.create_group('Space')
            spaceParameters.create_dataset('size', data=self.space.size)
            spaceParameters.create_dataset('material', data=self.materials.name(self.space.material))
            for subject in self.space.subjects:
                subjectParameters = spaceParameters.create_group(subject.__class__.__name__)
                for parameter_name, value in subject.__dict__.items():
                    subjectParameters.create_dataset(parameter_name, data=value)
            sourceParameters = group.create_group('Source')
            for parameter_name, value in self.source.__dict__.items():
                if parameter_name not in ('particles_emitted', 'timer', 'coordinates_table', 'size', 'rng_dist', 'rng_time', 'rng_dir', 'rng_ddist'):
                    sourceParameters.create_dataset(parameter_name, data=value)
            if self.subject is not None:
                subject = self.subject.__class__.__name__
                group.create_dataset('Subject', data=subject)
            group.create_dataset('Processes', data=Photons.processes)
            file.close()


class ParticleFlow(Process):
    """ Класс потока частиц """

    def __init__(self, source, space, materials, stop_time, particles_number, queue, solid_angle=None):
        super().__init__()
        self.source = source
        self.space = space
        self.materials = materials
        self.stop_time = stop_time
        self.particles_number = particles_number
        self.solid_angle = solid_angle
        self.queue = queue
        self.step = 1
        self.min_energy = 0
        self.daemon = True

    def off_the_solid_angle(self):
        if self.solid_angle is None:
            return []
        vector, angle = self.solid_angle
        cos_alpha = vector[0]*self.particles.direction[:, 0]
        cos_alpha += vector[1]*self.particles.direction[:, 1]
        cos_alpha += vector[2]*self.particles.direction[:, 2]
        indices = (cos_alpha <= cos(angle)).nonzero()[0]
        self.particles.delete(indices)
        return indices

    def invalid_particles(self, particles):
        indices = []
        indices.append((particles.energy <= self.min_energy).nonzero()[0])
        indices.append(self.space.outside(particles.coordinates))
        indices = np.concatenate(indices)
        return indices

    def valid_particles(self, particles):
        indices = self.space.inside(particles.coordinates)
        return indices

    def next_step(self, particles, interaction):
        free_path = interaction.get_free_path()
        particles.move(free_path)
        data = {}
        interaction_data = interaction.apply(self.valid_particles(particles))
        data.update({'Interaction': interaction_data})
        data.update({'Source timer': self.source.timer})
        data.update({'Step': self.step})
        self.queue.put(data)
        invalid_particles = self.invalid_particles(particles)
        invalid_particles_number = invalid_particles.size
        if self.source.timer <= self.stop_time:
            new_particles = self.source.generate_particles(invalid_particles_number)
            particles.replace(new_particles, invalid_particles)
        else:
            particles.delete(invalid_particles)
        self.step += 1

    def run(self):
        """ Реализация работы процесса """
        print(f'Flow started')
        particles = self.source.generate_particles(self.particles_number)
        interaction = Interaction(particles, self.space, self.materials)
        self.off_the_solid_angle()
        start = time()
        while particles.count > 0:
                self.next_step(particles, interaction)
        self.queue.put('Finish')
        print(f'Flow finished\n'
            + f'\tTime passed: {time() - start} seconds')


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
        self.timer = 0.
        self.coordinates_table = self._generate_coordinates_table()
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

    def _generate_coordinates_table(self):
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

    def set_timer(self, timer):
        self.timer = float(timer)

    def generate_coordinates(self, n):
        p = self.distribution.ravel()
        coordinates = self.coordinates_table[self.rng_dist.choice(np.arange(p.size), n, p=p)]
        dcoordinates = self.rng_ddist.uniform(0, self.voxel_size, coordinates.shape)
        coordinates += dcoordinates
        if self.rotated:
            coordinates -= self.rotation_center
            np.dot(coordinates, np.transpose(self.R), out=coordinates)
            coordinates += self.rotation_center
        coordinates += self.coordinates
        return coordinates

    def generate_emission_time(self, n):
        dt = log((self.nuclei_number + n)/self.nuclei_number)*self.half_life/log(2)
        a = 2**(-self.timer/self.half_life)
        b = 2**(-(self.timer + dt)/self.half_life)
        alpha = (b - a)*self.rng_time.random(n) + a
        emission_time = -log(alpha)*self.half_life/log(2)
        return emission_time, dt

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
        directions = self.generate_directions(n)
        coordinates = self.generate_coordinates(n)
        emission_time, dt = self.generate_emission_time(n)
        self.timer += dt
        particles = Photons(energies, directions, coordinates, emission_time)
        return particles

