import numpy as np
from subjects import Space
from phantoms import ae3cut
from collimators import SiemensSymbiaTSeriesLEHR
from detectors import SiemensSymbiaTSeries3_8
from modeling import Modeling
from modelingManagers import SourceManager
from materials import Materials
from hepunits import*


def main():
    angles = np.linspace(-pi/4, 3*pi/4, 32)
    projection_time = 15.*s
    pause_time = 1.*s

    materials = {
        'Compounds and mixtures/Air, Dry (near sea level)':         0,
        'Compounds and mixtures/Lung':                              1,
        'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
        'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
        'Compounds and mixtures/Sodium Iodide':                     4,
        'Elemental media/Pb':                                       5,
    }

    space = Space(
        size=(51.2*cm, 40.*cm, 60.*cm),
        material=0
    )

    detector = SiemensSymbiaTSeries3_8(
        coordinates=(0.*cm, 0.*cm, 0.*cm),
        size=space.size[:2]
    )

    collimator = SiemensSymbiaTSeriesLEHR(
        coordinates=(detector.coordinates[0], detector.coordinates[1], detector.size[2] + 0.5*cm),
        size=detector.size[:2]
    )

    phantom = ae3cut(
        coordinates=(collimator.coordinates[0], collimator.coordinates[1], collimator.coordinates[1] + collimator.size[2]),
    )

    source = SourceManager().efg3cut(
        coordinates=phantom.coordinates,
        activity=300*MBq,
    )

    space.add_subject(phantom)  
    space.add_subject(collimator)  
    space.add_subject(detector)

    materials = Materials(materials, max_energy=source.energy)

    materials.table = np.array([7, 7, 7, 10, 32, 82])

    for angle in angles:
        phantom.rotate((0., angle, 0.))
        source.rotate((0., angle, 0.))
        
        modeling = Modeling(
            space,
            source,
            materials,
            stop_time=source.timer + projection_time,
            particles_number=10**8,
            flow_number=10,
            file_name=f'efg3cut/{round(angle/deg, 1)} deg.hdf',
            iteraction_buffer=10**4,
            subject=detector
            )
        
        modeling.start()
        modeling.join()
        source.set_state(source.timer + pause_time)


if __name__ == '__main__':
    main()

