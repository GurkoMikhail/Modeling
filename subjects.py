import numpy as np
from numpy import sqrt, abs, mod, uint8, uint64
from utilites import culculate_R_euler, rotate_coordinates


class Space:
    """ Класс пространства моделирования """

    def __init__(self, size, material, subjects=[], **kwds):
        """ Конструктор пространства """
        self.size = np.asarray(size)    #cm
        self.material = material
        self.subjects = subjects
        self.args = ['']

        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])

    def outside(self, coordinates):
        """ Список попавших внутрь пространства"""
        outside = (coordinates > self.size) + (coordinates < 0)
        indices = outside.any(axis=1).nonzero()[0]
        return indices

    def inside(self, coordinates):
        """ Список попавших внутрь пространства """
        inside = (coordinates <= self.size)*(coordinates >= 0)
        indices = inside.all(axis=1).nonzero()[0]
        return indices

    def get_material(self, coordinates):
        """ Получить список веществ """
        material = np.full(coordinates.shape[0], self.material, uint8)
        for subject in self.subjects:
            off_subjects = (material == self.material).nonzero()[0]
            coordinates_off_subjects = coordinates[off_subjects]
            inside_subject = subject.inside(coordinates_off_subjects)
            subjects_material = subject.get_material_indices(coordinates_off_subjects[inside_subject])
            material[off_subjects[inside_subject]] = subjects_material
        return material

    @property
    def materials_list(self):
        materials = []
        for subject in self.subjects:
            materials.extend(np.unique(subject.material))
        materials = np.asarray(materials)
        return np.unique(materials)

    def add_subject(self, subject):
        """ Добавить объект """
        self.subjects.append(subject)


class Subject:
    """
    Класс обьекта

    [coordinates = (x, y, z)] = cm

    [size = (dx, dy, dz)] = cm
    
    [euler_angles = (alpha, beta, gamma)] = radian
    
    [material] = uint[:]
    """

    def __init__(self, coordinates, size, material, euler_angles=None, rotation_center=None):
        self.coordinates = np.asarray(coordinates)
        self.size = np.asarray(size)
        self.material = material
        self.rotated = False
        if euler_angles is not None:
            self.rotate(euler_angles, rotation_center)

    def rotate(self, euler_angles, rotation_center=None):
        self.rotated = True
        self.euler_angles = np.asarray(euler_angles)
        if rotation_center is None:
            rotation_center = np.asarray(self.size/2)
        self.rotation_center = rotation_center
        self.R = np.asarray(culculate_R_euler(self.euler_angles))

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates -= self.coordinates
        if self.rotated:
            coordinates -= self.rotation_center
            rotate_coordinates(coordinates, self.R)
            # np.dot(coordinates, np.transpose(self.R), out=coordinates)
            coordinates += self.rotation_center

    def inside(self, coordinates):
        """ Список попавших внутрь объекта с преобразованием координат """
        self.convert_to_local_coordinates(coordinates)
        inside = (coordinates <= self.size)*(coordinates >= 0)
        indices = inside.all(axis=1).nonzero()[0]
        return indices

    def get_material_indices(self, coordinates):
        """ Получить индексы материала """
        material = np.full(coordinates.shape[0], self.material, dtype=uint8)
        return material


class Phantom(Subject):
    """
    Класс фантома
    
    [coordinates = (x, y, z)] = cm
    
    [material] = uint[:,:,:]
    
    [voxel_size] = cm
    """

    def __init__(self, coordinates, material, voxel_size, euler_angles=None, rotation_center=None):
        size = np.asarray(material.shape)*voxel_size
        super().__init__(coordinates, size, material, euler_angles, rotation_center)
        self.voxel_size = voxel_size

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
    
    [material] = uint
    """

    def __init__(self, coordinates, size, material, hole_diameter, septa, space_material, euler_angles=None, rotation_center=None):
        super().__init__(coordinates, size, material, euler_angles, rotation_center)
        self.hole_diameter = hole_diameter
        self.septa = septa
        self.space_material = space_material
        y_period = sqrt(3)/2*self.hole_diameter + self.septa
        x_period = sqrt((2*y_period)**2 - y_period**2)
        self.period = np.stack((x_period, y_period))

    def get_collimated(self, coordinates):
        collimated = np.zeros(coordinates.shape[0], dtype=uint8)
        corners = np.stack((self.period, self.period/2))
        coordinates = mod(coordinates[:, :2], self.period)
        coordinates -= corners[1]
        corners -= corners[1]
        abs(coordinates, out=coordinates)
        a = sqrt(3)/4
        for corner in corners:
            dcoordinates = coordinates - corner
            dcoordinates /= self.hole_diameter
            abs(dcoordinates, out=dcoordinates)
            dx = dcoordinates[:, 0]
            dy = dcoordinates[:, 1]
            collimated += (dy <= a)*(a*dx + dy/4 <= a/2)
        return collimated.nonzero()[0]

    def get_material_indices(self, coordinates):
        material = super().get_material_indices(coordinates)
        collimated = self.get_collimated(coordinates)
        material[collimated] = self.space_material
        return material


class Detector(Subject):
    """ Класс детектора """


subjects_list = {
    'Subject': Subject,
    'Phantom': Phantom,
    'Collimator': Collimator,
    'Detector': Detector
    }

