import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling

size = np.asarray((51.2, 58.2, 40.))
space = Space(size, 0)

phantom = np.load('Phantoms/efg3cut_phantom.npy')
# phantom = np.rot90(phantom, k=2)
phantom = Phantom(
    coordinates=(0., 7., 0.),
    material=phantom,
    voxel_size=0.4)
space.add_subject(phantom)

collimator = Collimator((0., 7., 0.), (51.2, 40., 3.5), (0, np.pi/2, 0), 5, 0.15, 0.02)
space.add_subject(collimator)

detector = Detector((0., 3., 0.), (51.2, 40., 3), (0, np.pi/2, 0), 4)
space.add_subject(detector)

# source = np.load('Phantoms/efg3cut_source.npy')
# source = np.rot90(source, k=2)
# activity = np.sum(source)
# activity *= (0.4**3)*(10**6)

source = np.ones((128, 128, 100))
activity = 10**8

source = Source(
    coordinates=(0., 7., 0.),
    space=space,
    activity=activity,
    distribution=source,
    voxel_size=0.4,
    radiation_type='Gamma',
    energy=140.*10**3,
    half_life=6*60*60
    )

modeling = Modeling(space, source)
modeling.source.timer = 30.
modeling.start(60, 0.1)

