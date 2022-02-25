import numpy as np
from scipy.interpolate import interp1d
from h5py import File
from hepunits import*

class Materials:

    materials_table_path = 'tables/materials.hdf'
    lacs_table_path = 'tables/attenuationCoefficients.hdf'
    base_name = 'NIST'
    processes_names = {
        'PhotoelectricEffect': 'Photoelectric absorption',
        'ComptonScattering': 'Incoherent scattering',
        'CoherentScattering': 'Coherent scattering'
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
        self.materials_dict = {}
        for material, index in self.indices_dict.items():
            self.materials_dict.update({index: material})

    def select_atom(self, material):
        Z = self.table[material]
        return Z

    def name(self, material_index):
        return self.materials_dict[material_index]

    def construct_lac_funtions(self, processes):
        lac_funtions = {}
        lacsFile = File(self.lacs_table_path, 'r')
        lacsGroup = lacsFile[f'{self.base_name}/Linear attenuation coefficients']
        for process in processes:
            process_lac = {}
            for material, index in self.indices_dict.items():
                materialGroup = lacsGroup[material]
                energy = np.array(materialGroup['Energy'])
                max_energy_index = np.searchsorted(energy, self.max_energy, side='right') + 1
                energy = energy[:max_energy_index]
                lac = np.array(materialGroup[self.processes_names[process.name]][:max_energy_index])
                process_lac.update({index: interp1d(energy, lac/cm)})
            lac_funtions.update({process.name: process_lac})
        process_lac = {}
        for material, index in self.indices_dict.items():
            materialGroup = lacsGroup[material]
            energy = np.array(materialGroup['Energy'])
            max_energy_index = np.searchsorted(energy, self.max_energy, side='right') + 1
            energy = energy[:max_energy_index]
            total_lac = 0
            for process in processes:
                    total_lac += np.array(materialGroup[self.processes_names[process.name]][:max_energy_index])
            process_lac.update({index: interp1d(energy, total_lac/cm)})
        lac_funtions.update({'Total': process_lac})
        lacsFile.close()
        return lac_funtions

