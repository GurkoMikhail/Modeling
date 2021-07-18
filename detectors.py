from subjects import Detector


class SiemensSymbiaTSeries3_8(Detector):
    """
    Кристалл детектора 3/8 Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm

    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian
    """

    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 0.95
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
        dz = 1.59
        size = [*size, dz]
        material = 4
        super().__init__(coordinates, size, material, rotation_angles, rotation_center)

        