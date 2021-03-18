import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
from materials import Materials
import cProfile


if __name__ == '__main__':
    size = np.asarray((53.3, 60., 38.7))
    space = Space(size, 0)

    phantom = np.load('Phantoms/ae3_fix.npy')
    phantom = Phantom(
        coordinates=(1.05, (12.4 - 10.) + 5, -3.),
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


    materials = {
        'Compounds and mixtures/Air, Dry (near sea level)':         0,
        'Compounds and mixtures/Lung':                              1,
        'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
        'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
        'Compounds and mixtures/Sodium Iodide':                     4,
        'Elemental media/Pb':                                       5,
    }

    materials = Materials(materials, max_energy=140500)

    modeling = Modeling(
        space,
        source,
        materials,
        solid_angle=((0, -1, 0), 10*np.pi/180),
        time_step=0.01,
        file_name='efg3_fix 0.0 deg.hdf'
        )

    # modeling.start((0, 10.05))

    cProfile.run("modeling.start((0., 0.01))", 'efg3_fix 0.0 deg stats.txt')

