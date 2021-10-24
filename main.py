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
    angles = np.linspace(-np.pi/4, np.pi/2, 4)

    materials = {
        'Compounds and mixtures/Air, Dry (near sea level)':         0,
        'Compounds and mixtures/Lung':                              1,
        'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
        'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
        'Compounds and mixtures/Sodium Iodide':                     4,
        'Elemental media/Pb':                                       5,
    }

    space = Space(
        size=(53.3, 51.2, 60.),
        material=0
        )

    phantom = ae3(
        coordinates=((space.size[0] - 51.2)/2, (space.size[1] - 51.2)/2, 3.85),
        )
    space.add_subject(phantom)

    detector = SiemensSymbiaTSeries3_8(
        coordinates=((space.size[0] - 53.3)/2, (space.size[1] - 38.7)/2, 0),
        size=(53.3, 38.7)
        )
    space.add_subject(detector)

    collimator = SiemensSymbiaTSeriesLEHR(
        coordinates=(detector.coordinates[0], detector.coordinates[1], detector.size[2] + 0.5),
        size=detector.size[:2]
        )
    space.add_subject(collimator)

    source = efg3(
        coordinates=phantom.coordinates,
        activity=300*10**6,
        )

    materials = Materials(materials, max_energy=source.energy)

    materials.table = np.array([7, 7, 7, 10, 32, 82])

    for angle in angles:
        phantom.rotate((0., angle, 0.))
        source.rotate((0., angle, 0.))
        
        modeling = Modeling(
            space,
            source,
            materials,
            stop_time=0.1,
            particles_number=10**6,
            flow_number=2,
            file_name=f'efg3_full_angle {round(angle*180/np.pi, 1)} deg.hdf',
            subject=detector
            )
        
        modelings.append(modeling)
        modeling.start()

    for modeling in modelings:
        modeling.join()

