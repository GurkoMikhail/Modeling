import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
import multiprocessing as mp

def start_new_projection(angles, time):
    size = np.asarray((53.3, 70., 38.7))
    space = Space(size, 0)

    phantom = np.load('Phantoms/ae3_fix.npy')
    phantom = Phantom(
        coordinates=(1.05, (12.4 - 10.) + 10, -3.),
        material=phantom,
        voxel_size=0.4
        )
    space.add_subject(phantom)

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
        euler_angles=detector.euler_angles,
        rotation_center=detector.rotation_center
        )
    space.add_subject(collimator)

    source = Source(
        coordinates=phantom.coordinates,
        activity=300*10**6,
        distribution=np.load('Phantoms/efg3_fix.npy'),
        voxel_size=0.4,
        radiation_type='Gamma',
        energy=140.5*10**3,
        half_life=6*60*60
        )

    modeling = Modeling(
        space,
        source,
        solid_angle=((0, -1, 0), 15*np.pi/180),
        time_step=0.01,
        subject=detector
        )
    
    while not angles.empty():
        angle = angles.get()
        phantom.rotate((angle, 0, 0))
        source.rotate((angle, 0, 0))
        modeling.file_name = f'efg3_fix {round(angle*180/np.pi, 1)} deg.hdf'
        modeling.start(time)

if __name__ == '__main__':
    time = 15.
    angles = np.linspace(np.pi/4, -3*np.pi/4, 32)
    processes_number = 32

    queue = mp.Queue()
    for angle in angles:
        queue.put(angle)
    angles = queue

    processes = []
    for i in range(processes_number):
        process = mp.Process(target=start_new_projection, args=(angles, time))
        processes.append(process)
        process.start()

    for process in processes:
        process.join()
    print('End!')