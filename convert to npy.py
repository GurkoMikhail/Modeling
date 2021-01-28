import numpy as np

size = np.array([
    128,
    128,
    128
])

data = []
for x in open("Dat phantoms/ae3.dat"):
    x = float(x)
    # data.append(x)
    if x == 0.04:
        data.append(1)
    elif x == 0.15:
        data.append(2)
    elif x == 0.28:
        data.append(3)
    else:
        data.append(0)
data = np.array(data, dtype=np.uint8).reshape(size)
data = np.rot90(data, axes=(0, 2))
data = np.rot90(data, axes=(0, 1))
print(np.unique(data))
np.save('Phantoms/ae3.npy', data)