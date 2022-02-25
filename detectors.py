from subjects import Detector
from hepunits import*


class SiemensSymbiaTSeries3_8(Detector):
    """
    Кристалл детектора 3/8 Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm

    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian
    """

    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 9.5*mm
        size = [*size, dz]
        material = 4
        super().__init__(coordinates, size, material, rotation_angles, rotation_center)

        
class SiemensSymbiaTSeries5_8(Detector):
    """
    Кристалл детектора 5/8 Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm

    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian
    """

    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 15.9*mm
        size = [*size, dz]
        material = 4
        super().__init__(coordinates, size, material, rotation_angles, rotation_center)

        