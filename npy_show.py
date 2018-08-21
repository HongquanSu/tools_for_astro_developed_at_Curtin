import numpy as np
import sys

# have a quick look of a file containing numpy array that you saved
file = sys.argv[1]
array = np.load(file)
print(array)
print('shape:', np.shape(array))
