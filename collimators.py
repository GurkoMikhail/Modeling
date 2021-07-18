from subjects import Collimator

class SiemensSymbiaTSeriesLEHR(Collimator):
    """
    LEHR коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 2.405
        size = [*size, dz]
        material = 5
        hole_diameter = 0.111
        septa = 0.016
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        
class SiemensSymbiaTSeriesLEAP(Collimator):
    """
    LEAP коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 2.405
        size = [*size, dz]
        material = 5
        hole_diameter = 0.145
        septa = 0.02
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        
class SiemensSymbiaTSeriesLEUHR(Collimator):
    """
    LEUHR коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 3.58
        size = [*size, dz]
        material = 5
        hole_diameter = 0.116
        septa = 0.013
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        
class SiemensSymbiaTSeriesME(Collimator):
    """
    ME коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 4.064
        size = [*size, dz]
        material = 5
        hole_diameter = 0.294
        septa = 0.114
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        
class SiemensSymbiaTSeriesHE(Collimator):
    """
    HE коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 5.97
        size = [*size, dz]
        material = 5
        hole_diameter = 0.4
        septa = 0.2
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        