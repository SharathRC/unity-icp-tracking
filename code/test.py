import numpy as np

x = np.array([1,2,3,4,5,6,7,8,9]).reshape(3,3)
print(x)

print(np.mean(x, axis=0))