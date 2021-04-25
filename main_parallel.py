import os
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'

import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
from materials import Materials
import multiprocessing as mp

def start_new_projection(angles, time):
    size = np.asarray((51.2, 70., 40.))
    space = Space(size, 0)

    phantom = np.load('Phantoms/ae3cut.npy')
    phantom = Phantom(
        coordinates=(0., 12.4, 0.),
        material=phantom,
        voxel_size=0.4
        )
    space.add_subject(phantom)

    detector = Detector(
        coordinates=(0., 9.5, 0.),
        size=(51.2, 40., 9.5),
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
        coordinates=phantom.coordinates,
        activity=300*10**6,
        distribution=np.load('Phantoms/efg3cut.npy'),
        voxel_size=0.4,
        radiation_type='Gamma',
        energy=140.5*10**3,
        half_life=6*60*60
        )

    materials = {
        'Compounds and mixtures/Air, Dry (near sea level)':         0,
        'Compounds and mixtures/Lung':                              1,
        'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
        'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
        'Compounds and mixtures/Sodium Iodide':                     4,
        'Elemental media/Pb':                                       5,
    }

    materials = Materials(materials, max_energy=140500)

    import particles
    particles.Photons.processes.append('CoherentScattering')
    materials.table = np.array([7, 7, 7, 10, 32, 82])

    modeling = Modeling(
        space,
        source,
        materials,
        solid_angle=((0, -1, 0), 20*np.pi/180),
        time_step=0.01,
        subject=detector
        )
    
    while not angles.empty():
        angle = angles.get()
        phantom.rotate((angle, 0, 0))
        source.rotate((angle, 0, 0))
        modeling.file_name = f'efg3cut {round(angle*180/np.pi, 1)} deg.hdf'
        modeling.startMP(0., time, 9)

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