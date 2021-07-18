from modeling import Source
from numpy import load


class PointSource(Source):
    """
    Источник 99mТс-MIBI

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    
    [energy] = eV
    """

    def __init__(self, coordinates, activity, energy):
        distribution = [[[1.]]]
        voxel_size = 0.1
        radiation_type = 'Gamma'
        half_life = 6*60*60
        rotation_angles = None
        rotation_center = None
        super().__init__(coordinates, activity, distribution, voxel_size, radiation_type, energy, half_life, rotation_angles, rotation_center)

class Тс99m_MIBI(Source):
    """
    Источник 99mТс-MIBI

    [coordinates = (x, y, z)] = cm

    [activity] = Bq

    [distribution] = float[:,:,:]

    [voxel_size] = cm
    """

    def __init__(self, coordinates, activity, distribution, voxel_size, rotation_angles=None, rotation_center=None):
        radiation_type = 'Gamma'
        energy = 140.5*10**3
        half_life = 6*60*60
        super().__init__(coordinates, activity, distribution, voxel_size, radiation_type, energy, half_life, rotation_angles, rotation_center)


class efg3(Тс99m_MIBI):
    """
    Источник efg3

    [coordinates = (x, y, z)] = cm

    [activity] = Bq
    """

    def __init__(self, coordinates, activity, rotation_angles=None, rotation_center=None):
        distribution = load('Phantoms/efg3.npy')
        voxel_size = 0.4
        super().__init__(coordinates, activity, distribution, voxel_size, rotation_angles, rotation_center)

