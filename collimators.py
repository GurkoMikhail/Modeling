from subjects import Collimator
from hepunits import*

class SiemensSymbiaTSeriesLEHR(Collimator):
    """
    LEHR коллиматор Siemens Symbia T Series

    [coordinates = (x, y, z)] = cm
    
    [size = (dx, dy)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian

    [rotation_center = (x, y, z)] = cm
    """
    def __init__(self, coordinates, size, rotation_angles=None, rotation_center=None):
        dz = 24.05*mm
        size = [*size, dz]
        material = 5
        hole_diameter = 1.11*mm
        septa = 0.16*mm
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
        dz = 24.05*mm
        size = [*size, dz]
        material = 5
        hole_diameter = 1.45*mm
        septa = 0.2*mm
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
        dz = 35.8*mm
        size = [*size, dz]
        material = 5
        hole_diameter = 1.16*mm
        septa = 0.13*mm
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
        dz = 40.64*mm
        size = [*size, dz]
        material = 5
        hole_diameter = 2.94*mm
        septa = 1.14*mm
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
        dz = 59.7*mm
        size = [*size, dz]
        material = 5
        hole_diameter = 4.*mm
        septa = 2.*mm
        space_material = 0
        super().__init__(coordinates, size, material, hole_diameter, septa, space_material, rotation_angles, rotation_center)

        