from h5py._hl import group
import numpy as np
from h5py import File

class Materials:

    materials_table_path = 'tables/materials.hdf'
    lacs_table_path = 'tables/attenuationCoefficients.hdf'
    base_name = 'NIST'
    processes_names = {
        'PhotoelectricEffect': 'Photoelectric absorption',
        'ComptonScattering': 'Incoherent scattering'
    }

    def __init__(self, indices_dict, **kwds):
        self.indices_dict = indices_dict
        self.max_energy = 10**11
        self.args = [
            'max_energy'
            ]
        for arg in self.args:
            if arg in kwds:
                setattr(self, arg, kwds[arg])
        self._construct_tables()
        self.materials_dict = {}
        for material, index in self.indices_dict.items():
            self.materials_dict.update({index: material})

    def name(self, material_index):
        return self.materials_dict[material_index]

    def _construct_tables(self):
        self.lac_tables = [None]*len(self.indices_dict)
        lacsFile = File(self.lacs_table_path, 'r')
        lacsGroup = lacsFile[f'{self.base_name}/Linear attenuation coefficients']
        for material, index in self.indices_dict.items():
            lac = {}
            materialGroup = lacsGroup[material]
            energy = materialGroup['Energy']
            max_energy_index = np.searchsorted(energy, self.max_energy/10**6, side='right')
            for attribute, array in materialGroup.items():
                lac.update({attribute: array[:max_energy_index]})
            lac['Energy'] = lac['Energy']*10**6
            self.lac_tables[index] = lac
        lacsFile.close()

    def get_lac(self, materials, energy, processes):
        lac_out = np.zeros((len(processes), energy.size))
        for material in np.unique(materials):
            lac = self.lac_tables[material]
            indices = np.nonzero(materials == material)[0]
            for i, process in enumerate(processes):
                processName = self.processes_names[process.__class__.__name__]
                lac_out[i, indices] = np.interp(energy[indices], lac['Energy'], lac[processName])
        return lac_out

    def get_max_lac(self, materials, energy, processes):
        total_lac = np.zeros((materials.size, energy.size))
        for i, material in enumerate(materials):
            lac = self.lac_tables[material]
            for process in processes:
                processName = self.processes_names[process.__class__.__name__]
                total_lac[i] += np.interp(energy, lac['Energy'], lac[processName])
        return np.sum(total_lac, axis=0)

# mac_table = list([
#     np.load('macs/Air, Dry (near sea level).npy'),      # 0 - Воздух
#     np.load('macs/Lung Tissue.npy'),                    # 1 - Лёгочная ткань
#     np.load('macs/Tissue, Soft.npy'),                   # 2 - Мягкие ткани
#     np.load('macs/B-100 Bone-Equivalent Plastic.npy'),  # 3 - Эквивалент кости
#     np.load('macs/Sodium Iodide.npy'),                  # 4 - Йодит натрия
#     np.load('macs/Lead.npy'),                           # 5 - Свинец
# ])


# density_table = np.array([
#     1.205E-03,      # 0 - Воздух
#     1.205E-03,      # 1 - Лёгочная ткань
#     1.060E+00,      # 2 - Мягкие ткани
#     1.450E+00,      # 3 - Эквивалент кости
#     3.667E+00,      # 4 - Йодит натрия
#     1.135E+01,      # 5 - Свинец
# ])


# lac_table = list()
# for i, mac in enumerate(mac_table):
#     lac = mac
#     lac[:, 0] *= 10**6
#     lac[:, 1:] *= density_table[i]
#     lac_table.append(lac)


# materials_list = {
#     'Air, Dry (near sea level)':        0,
#     'Lung Tissue':                      1,
#     'Tissue, Soft':                     2,
#     'B-100 Bone-Equivalent Plastic':    3,
#     'Sodium Iodide':                    4,
#     'Lead':                             5,
# }

# materials_reverse_list = {
#     0:        'Air, Dry (near sea level)',
#     1:                      'Lung Tissue',
#     2:                     'Tissue, Soft',
#     3:    'B-100 Bone-Equivalent Plastic',
#     4:                    'Sodium Iodide',
#     5:                             'Lead',
# }


# def get_lac(materials, energy, processes):
#     lac_out = np.zeros((len(processes), energy.size))
#     for material in np.unique(materials):
#         lac = lac_table[material]
#         indices = np.nonzero(materials == material)[0]
#         for i, process in enumerate(processes):
#             process_index = process_indices[process.__class__.__name__]
#             lac_out[i, indices] = np.interp(energy[indices], lac[:, 0], lac[:, process_index])
#     return lac_out

# def get_max_lac(materials, energy, processes):
#     total_lac = np.zeros((materials.size, energy.size))
#     for i, material in enumerate(materials):
#         lac = lac_table[material]
#         for process in processes:
#             process_index = process_indices[process.__class__.__name__]
#             total_lac[i] += np.interp(energy, lac[:, 0], lac[:, process_index])
#     return np.sum(total_lac, axis=0)

