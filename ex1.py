import numpy as np
from numpy import linalg as la
M = ((1,-3,3), (3,-5,3), (6,-6,4))
T = np.zeros(3,3)
#T =M
print (T)
vals, vects = la.eig(M)
maxcol = list(vals).index(max(vals))
eigenvect = vects[:,maxcol]
print (eigenvect)