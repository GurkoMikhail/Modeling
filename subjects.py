import numpy as np
from numpy import cos, sin, sqrt, abs, mod, uint8, uint64, inf


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

    def ray_casting(self, coordinates, direction):
        """ Алгоритм бросания лучей для определения объекта местонахождения и длины пути до столкновения """
        path_length = np.full((coordinates.shape[0], ), inf)
        current_subject = np.zeros_like(path_length, dtype=uint8)
        for subject_index, subject in enumerate(self.subjects, 1):
            distance, inside_subject = subject.path_casting(coordinates, direction)
            current_subject[inside_subject] = subject_index
            intersectional = (path_length > distance).nonzero()[0]
            path_length[intersectional] = distance[intersectional]
        return current_subject, path_length

    def get_material_of_subject(self, subject_index):
        material = np.zeros_like(subject_index)
        complex_subject = np.zeros_like(subject_index, dtype=bool)
        for index, subject in enumerate(self.subjects, 1):
            indices = (subject_index == index).nonzero()[0]
            material[indices], complex_subject[indices] = subject.heaviest_material
        complex_subject = complex_subject.nonzero()[0]
        return material, complex_subject

    def outside(self, coordinates):
        """ Список непопавших внутрь пространства"""
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
            coordinates_off_subjects = subject.convert_to_local_coordinates(coordinates[off_subjects])
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
    
    [rotation_angles = (alpha, beta, gamma)] = radian
    
    [material] = uint[:]
    """

    def __init__(self, coordinates, size, material, rotation_angles=None, rotation_center=None):
        self.coordinates = np.asarray(coordinates)
        self.size = np.asarray(size)
        self.material = material
        self.rotated = False
        if rotation_angles is not None:
            self.rotate(rotation_angles, rotation_center)
        self._culculate_equation_coeffieients()
        self._find_heaviest_material()

    def rotate(self, rotation_angles, rotation_center=None):
        self.rotated = True
        self.rotation_angles = np.asarray(rotation_angles)
        if rotation_center is None:
            rotation_center = np.asarray(self.size/2)
        self.rotation_center = rotation_center
        alpha, beta, gamma = self.rotation_angles
        self.R = np.asarray([
            [cos(alpha)*cos(beta),  cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma),    cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma) ],
            [sin(alpha)*cos(beta),  sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),    sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma) ],
            [-sin(beta),            cos(beta)*sin(gamma),                                       cos(beta)*cos(gamma)                                    ]
        ])
        self.R = self.R.T

    def _find_heaviest_material(self):
        if type(self.material) is int:
            self.heaviest_material = (self.material, False)
        else:
            self.heaviest_material = (self.material.max(), True)

    def _culculate_equation_coeffieients(self):
        normals = np.asarray([
            [-1.,  0.,  0.],
            [ 0., -1.,  0.],
            [ 0.,  0., -1.],
            [ 1.,  0.,  0.],
            [ 0.,  1.,  0.],
            [ 0.,  0.,  1.],
        ])
        normals = normals.T
        D = np.asarray([
            0.,
            0.,
            0.,
            self.size[0],
            self.size[1],
            self.size[2]
        ])
        self.normals, self.D = normals, D

    def path_casting(self, coordinates, direction, local=False):
        if not local:
            coordinates = self.convert_to_local_coordinates(coordinates)
            direction = self.convert_to_local_direction(direction)
        inside = self.inside(coordinates)
        distance = -np.matmul(coordinates, self.normals)
        distance += self.D
        distance /= np.matmul(direction, self.normals)
        distance[distance <= 0] = inf
        distance = distance.min(axis=1) + 10**(-8)
        return distance, inside

    def convert_to_local_coordinates(self, coordinates):
        """ Преобразовать в локальные координаты """
        coordinates = coordinates.copy()
        coordinates -= self.coordinates
        if self.rotated:
            coordinates -= self.rotation_center
            np.matmul(coordinates, self.R, out=coordinates)
            coordinates += self.rotation_center
        return coordinates

    def convert_to_local_direction(self, direction):
        """ Преобразовать направление в локальное """
        direction = direction.copy()
        if self.rotated:
            np.matmul(direction, self.R, out=direction)
        return direction

    def inside(self, coordinates, local=True):
        """ Список попавших внутрь """
        if not local:
            coordinates = self.convert_to_local_coordinates(coordinates)
        inside = (coordinates <= self.size)*(coordinates >= 0)
        indices = inside.all(axis=1).nonzero()[0]
        return indices

    def get_material_indices(self, coordinates, local=True):
        """ Получить индексы материала """
        if not local:
            coordinates = self.convert_to_local_coordinates(coordinates)
        material = np.full(coordinates.shape[0], self.material, dtype=uint8)
        return material


class Phantom(Subject):
    """
    Класс фантома
    
    [coordinates = (x, y, z)] = cm
    
    [material] = uint[:,:,:]
    
    [voxel_size] = cm
    """

    def __init__(self, coordinates, material, voxel_size, rotation_angles=None, rotation_center=None):
        size = np.asarray(material.shape)*voxel_size
        super().__init__(coordinates, size, material, rotation_angles, rotation_center)
        self.voxel_size = voxel_size

    def get_material_indices(self, coordinates, local=True):
        if not local:
            coordinates = self.convert_to_local_coordinates(coordinates)
        coordinates = coordinates/self.voxel_size
        coordinates = coordinates.astype(uint64, copy=False)
        material_indices = self.material[(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])]
        return material_indices


class Collimator(Subject):
    """
    Класс коллиматора

    [coordinates = (x, y, z)] = cm
    
    [rotation_angles = (alpha, beta, gamma)] = radian
    
    [size = (dx, dy, dz)] = cm
    
    [hole_diameter] = cm
    
    [septa] = cm
    
    [material] = uint
    """

    def __init__(self, coordinates, size, material, hole_diameter, septa, space_material, rotation_angles=None, rotation_center=None):
        super().__init__(coordinates, size, material, rotation_angles, rotation_center)
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

    def get_material_indices(self, coordinates, local=True):
        material = super().get_material_indices(coordinates, local)
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

