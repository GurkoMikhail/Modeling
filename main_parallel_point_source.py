import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
import multiprocessing as mp

def start_new_projection(distances, time):
    size = np.asarray((53.3, 80., 38.7))
    space = Space(size, 0)

    detector = Detector(
        coordinates=(0., 9.5, 0.),
        size=(53.3, 38.7, 9.5),
        material=4,
        euler_angles=(0, np.pi/2, 0),
        rotation_center=(0., 0., 0.)
        )
    space.add_subject(detector)

    collimator = Collimator(
        coordinates=(detector.coordinates[0], detector.coordinates[1] + 0.5 + 2.4, detector.coordinates[2]),
        size=(*detector.size[:2], 2.4),
        material=5,
        hole_diameter=0.111,
        septa=0.016,
        space_material=space.material,
        euler_angles=detector.euler_angles,
        rotation_center=detector.rotation_center
        )
    space.add_subject(collimator)

    source = Source(
        coordinates=(size[0]/2, collimator.coordinates[1], size[2]/2),
        activity=70*10**6,
        distribution=np.ones((1, 1, 1)),
        voxel_size=0.2,
        radiation_type='Gamma',
        energy=140.5*10**3,
        half_life=6*60*60
        )

    modeling = Modeling(
        space,
        source,
        solid_angle=((0, -1, 0), 20*np.pi/180),
        time_step=0.01,
        subject=detector
        )
    
    while not distances.empty():
        distance = distances.get()
        source.coordinates[1] += distance
        space.size[1] = source.coordinates[1] + 1
        modeling.file_name = f'Point source {round(distance, 1)} sm.hdf'
        modeling.start((0., time))

if __name__ == '__main__':
    time = 10.
    distances = np.linspace(2.5, 80., 32)
    processes_number = 32

    queue = mp.Queue()
    for distance in distances:
        queue.put(distance)
    distances = queue

    processes = []
    for i in range(processes_number):
        process = mp.Process(target=start_new_projection, args=(distances, time))
        processes.append(process)
        process.start()

    for process in processes:
        process.join()
    print('End!')