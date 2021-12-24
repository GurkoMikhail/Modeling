from h5py import File
import numpy as np
from numpy import pi, sqrt, cos, sin, log
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
        self.solid_angle = None
        self.particles_number = 10**6
        self.iteraction_buffer = 10**3
        self.flow_number = 1
        self.file_name = f'{self}'
        self.subject = None
        self.save_emission_data = True
        self.save_dose_data = True
        self.distribution_voxel_size = 0.4
        self.modeling_data = {
            'Interactions data': {
                'Coordinates': [],
                'Energy transfer': [],
                'Emission time': [],
                'Emission coordinates': [],
                'Distance traveled': []
            },
            'Source timer': 0,
            'Events number': 0
        }
        self.args = [
            'start_time',
            'solid_angle',
            'particles_number',
            'iteraction_buffer',
            'flow_number',
            'file_name',
            'subject',
            'save_emission_data',
            'save_dose_data',
            'distribution_voxel_size'
            ]

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def generate_flows(self):
        queue = Queue(maxsize=1)
        flows = []
        for _ in range(self.flow_number):
            flow = ParticleFlow(
                source=self.source,
                space=self.space,
                materials=self.materials,
                stop_time=self.stop_time,
                particles_number=self.particles_number//self.flow_number,
                queue=queue,
                solid_angle=self.solid_angle
                )
            flows.append(flow)
        return flows, queue

    def start(self):
        if self.source._popen is None:
            self.source.start()
        return super().start()

    def run(self):
        print(f'{self.name} started')
        source_state = self.check_progress_in_file()
        self.source.set_state(*source_state)
        self.save_modeling_parameters()
        flows, queue = self.generate_flows()
        for flow in flows:
            flow.start()
        start_time = time()
        finished_flows = 0
        for step_data in iter(queue.get, None):
            finish_time = time() - start_time
            if step_data == 'Finish':
                finished_flows += 1
                if finished_flows == self.flow_number:
                    self.save_modeling_data()
                    print(f'\tReal time passed: {finish_time} seconds')
                    print(f'{self.name} ended!')
                    break
                continue
            self.update_modeling_data(step_data)
            if self.modeling_data['Events number'] >= self.iteraction_buffer:
                self.save_modeling_data()
                print(f'\tReal time passed: {finish_time} seconds')
                start_time = time()

    def check_progress_in_file(self):
        try:
            file = File(f'Output data/{self.file_name}', 'r')
            last_time = file['Source timer']
            last_time = float(np.array(last_time))
            state = None
            file.close()
        except Exception:
            print(f'\tНе удалось проверить прогресс')
            last_time = 0
            state = None
        finally:
            print(f'\tSource timer: {last_time}')
            return last_time, state

    def update_modeling_data(self, step_data):
        data = self.modeling_data['Interactions data']
        for dat in step_data['Interaction'].values():
            coordinates = dat['Coordinates']
            if len(coordinates) == 0:
                continue
            if self.save_dose_data:
                dose_distribution = np.histogramdd(
                    sample=coordinates,
                    bins=(self.space.size/self.distribution_voxel_size).astype(int),
                    range=((0, self.space.size[0]), (0, self.space.size[1]), (0, self.space.size[2])),
                    weights=dat['Energy transfer']
                )[0]
                if 'Dose distribution' in self.modeling_data:
                    self.modeling_data['Dose distribution'] += dose_distribution
                else:
                    self.modeling_data.update({'Dose distribution': dose_distribution})
            if self.save_emission_data:
                emission_distribution = np.histogramdd(
                    sample=dat['Emission coordinates'],
                    bins=(self.space.size/self.distribution_voxel_size).astype(int),
                    range=((0, self.space.size[0]), (0, self.space.size[1]), (0, self.space.size[2]))
                )[0]
                if 'Emission distribution' in self.modeling_data:
                    self.modeling_data['Emission distribution'] += emission_distribution
                else:
                    self.modeling_data.update({'Emission distribution': emission_distribution})
            if self.subject is not None:
                coordinates = self.subject.convert_to_local_coordinates(coordinates)
                inside_subject = self.subject.inside(coordinates)
                dat['Coordinates'] = coordinates
                for key in data.keys():
                    data[key].append(dat[key][inside_subject])
            else:
                for key in data.keys():
                    data[key].append(dat[key])
            self.modeling_data['Events number'] += data['Energy transfer'][-1].size
        self.modeling_data['Source timer'] = max([self.modeling_data['Source timer'], step_data['Source timer']])

    def concatenate_modeling_data(self):
        data = self.modeling_data['Interactions data']
        for key in data.keys():
            data[key] = np.concatenate(data[key])

    def save_modeling_data(self):
        self.concatenate_modeling_data()
        source_timer = self.modeling_data['Source timer']
        data = self.modeling_data['Interactions data']
        events_number = self.modeling_data['Events number']
        try:
            file = File(f'Output data/{self.file_name}', 'r+')
        except Exception:
            print(f'Не удалось сохранить {self.name}!')
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
                timer[...] = max([source_timer, float(np.array(timer))])
                group = file['Interactions data']
                for key in group.keys():
                    group[key].resize(
                        (group[key].shape[0] + data[key].shape[0]),
                        axis=0
                        )
                    group[key][-data[key].shape[0]:] = data[key]
            if self.save_dose_data:
                if not 'Dose distribution' in file:
                    group = file.create_group('Dose distribution')
                    volume = group.create_dataset('Volume', data=np.zeros((self.space.size/self.distribution_voxel_size).astype(np.uint), dtype=np.float64))
                    group.create_dataset('Voxel size', data=self.distribution_voxel_size)
                else:
                    group = file['Dose distribution']
                    volume = group['Volume']
                volume[...] += self.modeling_data['Dose distribution']
            if self.save_emission_data:
                if not 'Emission distribution' in file:
                    group = file.create_group('Emission distribution')
                    volume = group.create_dataset('Volume', data=np.zeros((self.space.size/self.distribution_voxel_size).astype(np.uint), dtype=np.float64))
                    group.create_dataset('Voxel size', data=self.distribution_voxel_size)
                else:
                    group = file['Emission distribution']
                    volume = group['Volume']
                volume[...] += self.modeling_data['Emission distribution']
            self.modeling_data = {
                'Interactions data': {
                    'Coordinates': [],
                    'Energy transfer': [],
                    'Emission time': [],
                    'Emission coordinates': [],
                    'Distance traveled': []
                },
                'Source timer': 0,
                'Events number': 0
            }
            print(f'{self.name} generated {events_number} events\n'
                + f'\tSource timer: {source_timer} seconds'
                )

    def save_modeling_parameters(self):
        try:
            file = File(f'Output data/{self.file_name}', 'x')
            group = file.create_group('Modeling parameters')
        except Exception:
            print(f'\tНе удалось записать параметры {self.name}')
        else:
            if self.solid_angle is not None:
                solidAngle = group.create_group('Solid angle')
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
            for parameter_name in dir(self.source):
                if parameter_name in ('distribution', 'coordinates', 'energy', 'voxel_size', 'initial_activity', 'radiation_type', 'half_life'):
                    sourceParameters.create_dataset(parameter_name, data=getattr(self.source, parameter_name))
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

    def off_the_solid_angle(self, direction):
        if self.solid_angle is None:
            return np.array([], dtype=int)
        vector, angle = self.solid_angle
        cos_alpha = vector[0]*direction[:, 0]
        cos_alpha += vector[1]*direction[:, 1]
        cos_alpha += vector[2]*direction[:, 2]
        indices = (cos_alpha <= cos(angle)).nonzero()[0]
        return indices

    def invalid_particles(self, particles):
        indices = []
        indices.append((particles.energy <= self.min_energy).nonzero()[0])
        indices.append(self.space.outside(particles.coordinates))
        indices.append(self.off_the_solid_angle(particles.direction))
        indices = np.concatenate(indices)
        return indices

    def send_data(self, data):
        data = {
            'Interaction': data,
            'Source timer': self.source.timer,
            'Step': self.step
            }
        self.queue.put(data)

    def next_step(self, particles, interaction):
        interaction_data = interaction.processing()
        self.send_data(interaction_data)
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
        print(f'\t{self.name} started')
        particles = self.source.generate_particles(self.particles_number)
        interaction = Interaction(particles, self.space, self.materials)
        start = time()
        while particles.count > 0:
                self.next_step(particles, interaction)
        self.queue.put('Finish')
        print(f'\t{self.name} finished\n'
            + f'\t\tTime passed: {time() - start} seconds')

