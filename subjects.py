from numpy import asarray, full, nonzero, uint8, uint64
from utilites import rotating_the_coordinates, culculate_R
from collimators import get_centers, get_collimated


class Space:
    """ Класс пространства моделирования """

    def __init__(self, size, material_index, subjects=[], **kwds):
        """ Конструктор пространства """
        self.size = asarray(size)    #cm
        self.material_index = material_index
        self.subjects = subjects
        self.args = ['']

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def get_material(self, coordinates):
        """ Получить список веществ """
        material = full(coordinates.shape[0], self.material_index, uint8)
        for subject in self.subjects:
            off_subjects = nonzero(material == self.material_index)[0]
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
        self.coordinates = asarray(coordinates)
        self.size = asarray(size)
        self.euler_angles = asarray(euler_angles)
        self.material_index = material_index
        self.R = asarray(culculate_R(self.euler_angles))

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates[:, 0] -= self.coordinates[0]
        coordinates[:, 1] -= self.coordinates[1]
        coordinates[:, 2] -= self.coordinates[2]
        rotating_the_coordinates(coordinates, self.R)

    def inside(self, coordinates):
        """ Список попавших внутрь объекта с преобразованием координат """
        self.convert_to_local_coordinates(coordinates)
        in_x = (coordinates[:, 0] <= self.size[0])*(coordinates[:, 0] >= 0)
        in_y = (coordinates[:, 1] <= self.size[1])*(coordinates[:, 1] >= 0)
        in_z = (coordinates[:, 2] <= self.size[2])*(coordinates[:, 2] >= 0)
        indices = nonzero(in_x*in_y*in_z)[0]
        return indices

    def get_material_indices(self, coordinates):
        """ Получить индексы материала """
        material = full(coordinates.shape[0], self.material_index, dtype=uint8)
        return material


class Phantom(Subject):
    """
    Класс фантома
    
    [coordinates = (x, y, z)] = cm
    
    [material] = uint[:,:,:]
    
    [voxel_size] = cm
    """

    def __init__(self, coordinates, material, voxel_size):
        self.coordinates = asarray(coordinates)
        self.material = asarray(material)
        self.voxel_size = voxel_size
        self.size = asarray(self.material.shape)*self.voxel_size

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates[:, 0] -= self.coordinates[0]
        coordinates[:, 1] -= self.coordinates[1]
        coordinates[:, 2] -= self.coordinates[2]

    def get_material_indices(self, coordinates):
        coordinates = coordinates/self.voxel_size
        coordinates = coordinates.astype(uint64, copy=False)
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
        self.holes_centers = get_centers(self.size, self.hole_diameter, self.septa)

    def get_material_indices(self, coordinates):
        material = full(coordinates.shape[0], self.material_index, dtype=uint8)
        collimated = get_collimated(coordinates, self.holes_centers, self.hole_diameter)
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

