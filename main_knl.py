from multiprocessing import Lock, Queue, Process
import numpy as np
from subjects import Space, Phantom, Collimator, Detector
from modeling import Source, Modeling
from materials import Materials

class ModelingProcess(Process):

    def __init__(self, parameters, lock, group=None, target=None, name=None, args=(), kwargs={}, *, daemon=None):
        super().__init__(group=group, target=target, name=name, args=args, kwargs=kwargs, daemon=daemon)
        self.parameters = parameters
        self.lock = lock

    def run(self):
        while not self.parameters.empty():
            parameters = self.parameters.get()
            detector = Detector(
                coordinates=(0., 9.5, 0.),
                size=(51.2, 40., 9.5),
                material=4,
                euler_angles=(0, np.pi/2, 0),
                rotation_center=(0., 0., 0.)
                )

            collimator = Collimator(
                coordinates=(detector.coordinates[0], detector.coordinates[1] + 0.5 + parameters['Collimator thickness'], detector.coordinates[2]),
                size=(*detector.size[:2], parameters['Collimator thickness']),
                material=5,
                hole_diameter=parameters['Hole diameter'],
                septa=parameters['Septa'],
                euler_angles=detector.euler_angles,
                rotation_center=detector.rotation_center
                )

            phantom = np.load('Phantoms/ae3cut.npy')
            phantom = Phantom(
                coordinates=(0., collimator.coordinates[1] - 5., 0.),
                material=phantom,
                voxel_size=0.4
                )

            size = np.asarray((detector.size[0], phantom.coordinates[1] + phantom.size[1], detector.size[1]))
            space = Space(size, 0)
            space.add_subject(detector)
            space.add_subject(collimator)
            space.add_subject(phantom)

            source = Source(
                coordinates=phantom.coordinates,
                activity=300*10**6,
                distribution=np.load('Phantoms/efg3cut.npy'),
                voxel_size=0.4,
                radiation_type='Gamma',
                energy=140.5*10**3,
                half_life=6*60*60
                )

            materials = {
                'Compounds and mixtures/Air, Dry (near sea level)':         0,
                'Compounds and mixtures/Lung':                              1,
                'Compounds and mixtures/Tissue, Soft (ICRU-44)':            2,
                'Compounds and mixtures/B-100 Bone-Equivalent Plastic':     3,
                'Compounds and mixtures/Sodium Iodide':                     4,
                'Elemental media/Pb':                                       5,
            }

            materials = Materials(materials, max_energy=140500)

            modeling = Modeling(
                space,
                source,
                materials,
                solid_angle=((0, -1, 0), 15*np.pi/180),
                time_step=0.01,
                subject=detector,
                file_name=parameters['Name'] + '.hdf'
                )
            
            modeling.startMP(parameters['Time'], self.lock)

def generateParameters():
    cores_number = 72
    totalTime = 15.
    startTime = np.round(np.linspace(0, totalTime, cores_number + 1), 2)[:-1]
    finishTime = np.round((startTime + (startTime[1] - startTime[0])), 2)
    times = np.stack((startTime, finishTime), axis=1)
    collimatorThickness = [2.405, 2.405, 3.58, 4.064]
    holeDiameter = [0.111, 0.145, 0.116, 0.294]
    septa = [0.016, 0.02, 0.013, 0.114]
    names = ['LEHR', 'LEAP', 'LEUHR', 'ME']
    parameters_stack = []
    for i, name in enumerate(names):
        for time in times:
            parameters = {
                'Time': time,
                'Collimator thickness': collimatorThickness[i],
                'Hole diameter': holeDiameter[i],
                'Septa': septa[i],
                'Name': name
                }
            parameters_stack.append(parameters)
    return parameters_stack

if __name__ == '__main__':
    process_number = 288
    parameters_stack = generateParameters()
    queue_stack = [Queue() for _ in range(process_number)]
    lock = Lock()
    processes_stack = []
    ratio = len(parameters_stack)//process_number
    for i, queue in enumerate(queue_stack):
        for parameters in parameters_stack[i*ratio:(i + 1)*ratio]:
            queue.put(parameters)
        processes_stack.append(ModelingProcess(queue, lock))
    for process in processes_stack:
        process.start()
    for process in processes_stack:
        process.join()
    print('End!')