import numpy as np
from subjects import Space, Collimator, Detector
from modeling import Source, Modeling
from materials import Materials


if __name__ == '__main__':
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
        voxel_size=0.001,
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
        solid_angle=((0, -1, 0), 45*np.pi/180),
        time_step=0.01,
        subject=detector
        )

    distance = 10

    source.coordinates[1] += distance
    space.size[1] = source.coordinates[1] + 1
    modeling.file_name = f'Point source {round(distance, 1)} sm.hdf'
    
    modeling.startMP(0, 0.1, 2)

