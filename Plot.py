
import numpy as np 
import math
from matplotlib import pyplot as plt
import copy
#from Fagverk import *

"""
fig = plt.figure()
x,y = nodes.T
plt.scatter(x,y)
plt.savefig('Konstruksjon.png')

#plt.show()

"""
a = np.array([1,2])
b = np.array([[0,-1], [1,0]])
c = np.matmul(a,b)
d = np.matmul(a,c)
print(d)