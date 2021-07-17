import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
from materials import Materials


if __name__ == '__main__':
    size = np.asarray((53.3, 60., 38.7))
    space = Space(size, 0)

    phantom = np.load('Phantoms/ae3.npy')
    phantom = Phantom(
        coordinates=(1.05, 12.4, -3.),
        material=phantom,
        voxel_size=0.4,
        # rotation_angles=(0, 0, np.pi/2)
        )
    space.add_subject(phantom)

    detector = Detector(
        coordinates=(0., 0.95, 0.),
        size=(53.3, 38.7, 0.95),
        material=4,
        rotation_angles=(0, 0, -np.pi/2),
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
        rotation_angles=detector.rotation_angles,
        rotation_center=detector.rotation_center
        )
    space.add_subject(collimator)

    source = Source(
        coordinates=phantom.coordinates,
        activity=300*10**6,
        distribution=np.load('Phantoms/efg3.npy'),
        voxel_size=0.4,
        radiation_type='Gamma',
        energy=140.5*10**3,
        half_life=6*60*60,
        # rotation_angles=phantom.rotation_angles,
        # rotation_center=phantom.rotation_center
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

    materials.table = np.array([7, 7, 7, 10, 32, 82])


    modeling = Modeling(
        space,
        source,
        materials,
        stop_time=1,
        particles_number=10**7,
        flow_number=2,
        file_name='efg3_full_angle 0.0 deg.hdf',
        subject=detector
        )

    modeling.start()
    modeling.join()

