from h5py import File
import numpy as np

if __name__ == '__main__':
    materialsFile = File('tables/materials.hdf', 'r')
    macFile = File('tables/mac.hdf', 'a')
    materialsGroup = materialsFile['NIST/Compounds and mixtures']
    elementalGroup = macFile['NIST/Elemental media']
    try:
        compoundsGroup = macFile['NIST/Compounds and mixtures']
    except Exception:
        compoundsGroup = macFile.create_group('NIST/Compounds and mixtures')
    for material, parameters in materialsGroup.items():
        macGroup = compoundsGroup.create_group(material)
        energy = None
        for element in parameters['Composition'].keys():
            newEnergy = elementalGroup[element]['Energy']
            if energy is None:
                energy = newEnergy
            mask = np.in1d(newEnergy, energy, invert=True)
            energy = np.concatenate((energy, newEnergy[mask]))
        energy = np.sort(energy)
        macGroup.create_dataset('Energy', data=energy)
        data = {
            'Coherent scattering': 0,
            'Incoherent scattering': 0,
            'Photoelectric absorption': 0,
            'Pair production in nuclear field': 0,
            'Pair production in electron field': 0
        }
        for element, value in parameters['Composition'].items():
            for effect, array in elementalGroup[element].items():
                if effect == 'Energy':
                    continue
                data[effect] += value*np.interp(energy, elementalGroup[element]['Energy'], array)
        for effect, array in data.items():
            macGroup.create_dataset(effect, data=array)
    materialsFile.close()
    macFile.close()

