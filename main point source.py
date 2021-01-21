import numpy as np
from subjects import Space, Subject, Phantom, Collimator, Detector
from modeling import Source

size = np.asarray((20., 15., 20.))
space = Space(size, 0)

collimator = Collimator((0., 7., 0.), (20., 20., 3.5), (0, np.pi/2, 0), 5, 0.15, 0.02)
space.add_subject(collimator)

detector = Detector((0., 3., 0.), (20., 10., 3), (0, np.pi/2, 0), 4)
space.add_subject(detector)

source = np.full((1, 1, 1), 1)
activity = 70*10**6

source = Source(
    coordinates=(size[0]/2, 7. + 5., size[2]/2),
    space=space,
    activity=activity,
    distribution=source,
    voxel_size=0.1,
    radiation_type='Gamma',
    energy=140.*10**3,
    half_life=6*60*60
    )

source.start(10, 0.1)

