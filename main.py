import numpy as np
from subjects import Space
from phantoms import ae3
from collimators import SiemensSymbiaTSeriesLEHR
from detectors import SiemensSymbiaTSeries3_8
from modeling import Modeling
from sources import efg3
from materials import Materials

if __name__ == '__main__':
    modelings = []
    angles = np.linspace(-np.pi/4, 3*np.pi/4, 60)

    size = np.asarray((53.3, 100., 53.3))
    space = Space(size, 0)

    phantom = ae3(
        coordinates=(1.05, 3.85, 1.05),
        )
    space.add_subject(phantom)

    detector = SiemensSymbiaTSeries3_8(
        coordinates=(0., 0.95, 0.),
        size=(53.3, 53.3),
        rotation_angles=(0, 0, -np.pi/2),
        rotation_center=(0., 0., 0.)
        )
    space.add_subject(detector)

    collimator = SiemensSymbiaTSeriesLEHR(
        coordinates=(detector.coordinates[0], detector.coordinates[1] + 0.5 + 2.4, detector.coordinates[2]),
        size=detector.size[:2],
        rotation_angles=detector.rotation_angles,
        rotation_center=detector.rotation_center
        )
    space.add_subject(collimator)

    source = efg3(
        coordinates=phantom.coordinates,
        activity=300*10**6,
        )

    materials = {
        'Compounds and mixtures/Air, Dry (near sea level)':         0,
        'Compounds and mixtures/Lung':                              1,
        'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
        'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
        'Compounds and mixtures/Sodium Iodide':                     4,
        'Elemental media/Pb':                                       5,
    }

    materials = Materials(materials, max_energy=source.energy)

    materials.table = np.array([7, 7, 7, 10, 32, 82])

    for angle in angles:
        phantom.rotate((0., 0., angle))
        source.rotate((0., 0., angle))
        
        modeling = Modeling(
            space,
            source,
            materials,
            stop_time=20.,
            particles_number=10**7,
            flow_number=1,
            file_name=f'efg3_full_angle {round(angle*180/np.pi, 1)} deg.hdf',
            subject=detector
            )
        
        modelings.append(modeling)
        modeling.start()

    for modeling in modelings:
        modeling.join()

