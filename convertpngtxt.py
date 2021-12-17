from pathlib import Path
import numpy as np
import cv2


directory = Path('X:/Caro')

for file in directory.glob('inputs/*.png'):
	try:
		image = cv2.imread(str(file.absolute()), cv2.IMREAD_GRAYSCALE)

		with open(file.absolute().parent / (file.stem + '.txt'), 'w') as fd:
			image = np.array(image)
			s = [str(image.shape[0]) + ' ' + str(image.shape[1])]
			for x in range(image.shape[0]):
				s.append(''.join([('0' if image[x,y]==0 else '1') for y in range(image.shape[1])]))
			fd.write('\n'.join(s))
	except Exception as e:
		print('Skipping file due to exception thrown: ' + file)
		print('Due to: ' + str(e))
		continue