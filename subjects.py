import numpy as np
import utilites
import collimators
from time import time

class Space:
    """ Класс пространства моделирования """

    def __init__(self, size, material_index, subjects=[], **kwds):
        """ Конструктор пространства """
        self.size = np.asarray(size)    #cm
        self.material_index = material_index
        self.subjects = subjects
        self.args = ['']

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def get_material(self, coordinates):
        """ Получить список веществ """
        material = np.full(coordinates.shape[0], self.material_index, np.uint8)
        for subject in self.subjects:
            off_subjects = np.nonzero(material == self.material_index)[0]
            coordinates_off_subjects = coordinates[off_subjects]
            inside_subject = subject.inside(coordinates_off_subjects)
            subjects_material = subject.get_material_indices(coordinates_off_subjects[inside_subject])
            material[off_subjects[inside_subject]] = subjects_material
        return material

    def add_subject(self, subject):
        """ Добавить объект """
        self.subjects.append(subject)


class Subject:
    """
    Класс обьекта

    [coordinates = (x, y, z)] = cm

    [size = (dx, dy, dz)] = cm
    
    [euler_angles = (alpha, beta, gamma)] = radian
    
    [material_index] = uint[:]
    """

    def __init__(self, coordinates, size, euler_angles, material_index):
        self.coordinates = np.asarray(coordinates)
        self.size = np.asarray(size)
        self.euler_angles = np.asarray(euler_angles)
        self.material_index = material_index
        self.R = np.asarray(utilites.culculate_R(self.euler_angles))

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates[:, 0] -= self.coordinates[0]
        coordinates[:, 1] -= self.coordinates[1]
        coordinates[:, 2] -= self.coordinates[2]
        utilites.rotating_the_coordinates(coordinates, self.R)

    def inside(self, coordinates):
        """ Список попавших внутрь объекта с преобразованием координат """
        self.convert_to_local_coordinates(coordinates)
        in_x = (coordinates[:, 0] <= self.size[0])*(coordinates[:, 0] >= 0)
        in_y = (coordinates[:, 1] <= self.size[1])*(coordinates[:, 1] >= 0)
        in_z = (coordinates[:, 2] <= self.size[2])*(coordinates[:, 2] >= 0)
        indices = np.nonzero(in_x*in_y*in_z)[0]
        return indices

    def get_material_indices(self, coordinates):
        """ Получить индексы материала """
        material = np.full(coordinates.shape[0], self.material_index, dtype=np.uint8)
        return material

class Phantom(Subject):
    """
    Класс фантома
    
    [coordinates = (x, y, z)] = cm
    
    [material] = uint[:,:,:]
    
    [voxel_size] = cm
    """

    def __init__(self, coordinates, material, voxel_size):
        self.coordinates = np.asarray(coordinates)
        self.material = np.asarray(material)
        self.voxel_size = voxel_size
        self.size = np.asarray(self.material.shape)*self.voxel_size

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates[:, 0] -= self.coordinates[0]
        coordinates[:, 1] -= self.coordinates[1]
        coordinates[:, 2] -= self.coordinates[2]

    def get_material_indices(self, coordinates):
        coordinates = coordinates/self.voxel_size
        coordinates = coordinates.astype(np.uint32, copy=False)
        material_indices = self.material[(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])]
        return material_indices


class Collimator(Subject):
    """
    Класс коллиматора

    [coordinates = (x, y, z)] = cm
    
    [euler_angles = (alpha, beta, gamma)] = radian
    
    [size = (dx, dy, dz)] = cm
    
    [hole_diameter] = cm
    
    [septa] = cm
    
    [material_index] = uint
    """

    def __init__(self, coordinates, size, euler_angles, material_index, hole_diameter, septa):
        super().__init__(coordinates, size, euler_angles, material_index)
        self.hole_diameter = hole_diameter
        self.septa = septa
        self.holes_centers = collimators.get_centers(self.size, self.hole_diameter, self.septa)

    def get_material_indices(self, coordinates):
        material = np.full(coordinates.shape[0], self.material_index, dtype=np.uint8)
        collimated = collimators.get_collimated(coordinates, self.holes_centers, self.hole_diameter)
        material[collimated] = 0
        return material


class Detector(Subject):
    """ Класс детектора """


subjects_list = {
    'Subject': Subject,
    'Phantom': Phantom,
    'Collimator': Collimator,
    'Detector': Detector
    }



# class GammaCamera(Collimator):
#     """ Класс гамма-камеры """

#     def __init__(self, coordinates, euler_angles, collimator_size, hole_diameter, septa, gap, detector_thickness, collimator_material_index, detector_material_index, **kwds):
#         self.coordinates = coordinates
#         self.euler_angles = np.asarray(euler_angles)
#         self.collimator_size = np.asarray(collimator_size)
#         self.hole_diameter = hole_diameter
#         self.septa = septa
#         self.gap = gap
#         self.detector_thickness = detector_thickness
#         self.collimator_material_index = collimator_material_index
#         self.detector_material_index = detector_material_index
#         self.size = np.asarray(collimator_size)
#         self.size[2] += gap + detector_thickness
#         self.holes_centers = self.get_centers()
#         self.R = utilites.culculate_R(self.euler_angles)
#         self.tops = utilites.find_tops(self.coordinates, self.size, self.euler_angles)
#         self.args = ['']

#         for arg in self.args:
#             if arg in kwds:
#                 setattr(self, arg, kwds[arg])    

#     def inside_detector(self, coordinates):
#         indices = np.nonzero(coordinates[:, 2] > self.collimator_size[2] + self.gap)[0]
#         return indices

#     def inside_collimator(self, coordinates):
#         indices = np.nonzero(coordinates[:, 2] <= self.collimator_size[2])
#         return indices

#     def get_material_indices(self, coordinates):
#         material = np.full(coordinates.shape[0], 1, dtype=np.uint8)
#         inside_detector = self.inside_detector(coordinates)
#         material[inside_detector] = self.detector_material_index
#         inside_collimator = self.inside_collimator(coordinates)
#         collimated = get_collimated(coordinates[inside_collimator], self.holes_centers, self.hole_diameter)
#         material[collimated] = self.collimator_material_index
#         return material

#     def get_data(self, coordinates):
#         inside_region = utilites.inside_region(coordinates, self.tops)
#         coordinates = self.convert_to_local_coordinates(coordinates[inside_region])
#         material = self.get_material_indices(coordinates)
#         return (inside_region, material)
