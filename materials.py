from numba import jit
from numba.typed import List
from numpy import array, load, unique, zeros, nonzero, interp

mac_table = List([
    load('macs/Air, Dry (near sea level).npy'),      # 0 - Воздух
    load('macs/Lung Tissue.npy'),                    # 1 - Лёгочная ткань
    load('macs/Tissue, Soft.npy'),                   # 2 - Мягкие ткани
    load('macs/B-100 Bone-Equivalent Plastic.npy'),  # 3 - Эквивалент кости
    load('macs/Sodium Iodide.npy'),                  # 4 - Йодит натрия
    load('macs/Lead.npy'),                           # 5 - Свинец
])

density_table = array([
    1.205E-03,      # 0 - Воздух
    1.205E-03,      # 1 - Лёгочная ткань
    1.060E+00,      # 2 - Мягкие ткани
    1.450E+00,      # 3 - Эквивалент кости
    3.667E+00,      # 4 - Йодит натрия
    1.135E+01,      # 5 - Свинец
])

lac_table = List()
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

process_indices = {
    'CoherentScatter': 1,
    'ComptonScattering': 2,
    'PhotoelectricEffect': 3
}

def get_lac(material, energy, processes):
    lac_out = zeros((len(processes), energy.size))
    for m in unique(material):
        lac = lac_table[m]
        indices = nonzero(material == m)[0]
        for i, process in enumerate(processes):
            process_index = process_indices[process.__class__.__name__]
            lac_out[i, indices] = interp(energy[indices], lac[:, 0], lac[:, process_index])
    return lac_out

        
