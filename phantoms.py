from subjects import Phantom
from numpy import load


class ae3(Phantom):
    """
    Карта ослабления ae3
    
    [coordinates = (x, y, z)] = cm
    """

    def __init__(self, coordinates, rotation_angles=None, rotation_center=None):
        material = load('Phantoms/ae3.npy')
        voxel_size = 0.4
        super().__init__(coordinates, material, voxel_size, rotation_angles, rotation_center)

