import matplotlib.pyplot as plt
import numpy as np
import cv2
from pathlib import Path


dbg_path = Path('./dbg')

base = cv2.imread(str(dbg_path / 'p0.png'), cv2.IMREAD_COLOR)
tranches = {int(p.stem.split('tranche')[1]): cv2.imread(str(p), cv2.IMREAD_COLOR) for p in dbg_path.glob('p1tranche*.png')}
print(list(tranches.keys()))

partition_size = 16

cum = np.zeros((partition_size, partition_size), dtype=int)
for tranche_number, tranche_img in tranches.items():
    dims = np.array(tranche_img.shape[:2]) // partition_size
    for outer_x in range(dims[0]):
        for outer_y in range(dims[1]):
            for x in range(partition_size):
                for y in range(partition_size):
                    if base[outer_x * partition_size + x, outer_y * partition_size + y, 0] != tranche_img[outer_x * partition_size + x, outer_y * partition_size + y, 0]:
                        cum[x,y] += 1
    print(f'Tranche number {tranche_number} was analyzed')

denom = np.zeros_like(cum)
dims = np.array(base.shape[:2]) // partition_size
for outer_x in range(dims[0]):
    for outer_y in range(dims[1]):
        for x in range(partition_size):
            for y in range(partition_size):
                if base[outer_x * partition_size + x, outer_y * partition_size + y, 0] != 0:
                    denom[x,y] += 1

cum = cum.astype(np.float32) / denom.astype(np.float32)
print(cum)
plt.imshow(cum)
plt.show()
