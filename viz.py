from matplotlib import pyplot as plt
from sys import stdin
import numpy as np


xs = []
ys = []
ignore = True
img = False
image = []
for line in stdin:
	if line.rstrip() == 'Path:':
		ignore = False
	elif line.rstrip() == 'Remaining:':
		img = True
	elif img:
		image.append([int(s) for s in line.split()])
	elif not ignore:
		line = line.split()
		xs.append(int(line[1])-1)
		ys.append(int(line[0])-1)
image = np.array(image)
print(image.shape)
#plt.plot(xs, ys)#, marker = '.')
plt.imshow(image)
plt.savefig('./prev_viz.png')
plt.show()