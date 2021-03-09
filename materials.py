import numpy as np


mac_table = list([
    np.load('macs/Air, Dry (near sea level).npy'),      # 0 - Воздух
    np.load('macs/Lung Tissue.npy'),                    # 1 - Лёгочная ткань
    np.load('macs/Tissue, Soft.npy'),                   # 2 - Мягкие ткани
    np.load('macs/B-100 Bone-Equivalent Plastic.npy'),  # 3 - Эквивалент кости
    np.load('macs/Sodium Iodide.npy'),                  # 4 - Йодит натрия
    np.load('macs/Lead.npy'),                           # 5 - Свинец
])


density_table = np.array([
    1.205E-03,      # 0 - Воздух
    1.205E-03,      # 1 - Лёгочная ткань
    1.060E+00,      # 2 - Мягкие ткани
    1.450E+00,      # 3 - Эквивалент кости
    3.667E+00,      # 4 - Йодит натрия
    1.135E+01,      # 5 - Свинец
])


lac_table = list()
for i, mac in enumerate(mac_table):
    lac = mac
    lac[:, 0] *= 10**6
    lac[:, 1:] *= density_table[i]
    lac_table.append(lac)


materials_list = {
    'Air, Dry (near sea level)':        0,
    'Lung Tissue':                      1,
    'Tissue, Soft':                     2,
    'B-100 Bone-Equivalent Plastic':    3,
    'Sodium Iodide':                    4,
    'Lead':                             5,
}

materials_reverse_list = {
    0:        'Air, Dry (near sea level)',
    1:                      'Lung Tissue',
    2:                     'Tissue, Soft',
    3:    'B-100 Bone-Equivalent Plastic',
    4:                    'Sodium Iodide',
    5:                             'Lead',
}

process_indices = {
    'CoherentScattering': 1,
    'ComptonScattering': 2,
    'PhotoelectricEffect': 3
}


def get_lac(materials, energy, processes):
    lac_out = np.zeros((len(processes), energy.size))
    for material in np.unique(materials):
        lac = lac_table[material]
        indices = np.nonzero(materials == material)[0]
        for i, process in enumerate(processes):
            process_index = process_indices[process.__class__.__name__]
            lac_out[i, indices] = np.interp(energy[indices], lac[:, 0], lac[:, process_index])
    return lac_out

def get_max_lac(materials, energy, processes):
    total_lac = np.zeros((materials.size, energy.size))
    for i, material in enumerate(materials):
        lac = lac_table[material]
        for process in processes:
            process_index = process_indices[process.__class__.__name__]
            total_lac[i] += np.interp(energy, lac[:, 0], lac[:, process_index])
    return np.sum(total_lac, axis=0)

