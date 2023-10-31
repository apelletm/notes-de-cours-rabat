from fpylll import *
import copy

from fpylll.algorithms.bkz import BKZReduction
from math import log
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
import matplotlib.pyplot as plt



"""
BKZ demo
"""


A = IntegerMatrix.random(80, "qary", k=40, bits=45)

LLL.reduction(A)
M = GSO.Mat(A)
M.update_gso()
d = A.nrows

colours = ["#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", "#F15854"]
norms = [[log(M.get_r(i,i)) for i in range(d)]]
plt.plot(norms[0],label="lll", color=colours[0])

beta = 40
tours = 8
par = BKZ.Param(block_size=beta,strategies=BKZ.DEFAULT_STRATEGY)
bkz = BKZ2(M)

for i in range(tours):
    bkz.tour(par)
    norms += [[log(M.get_r(j,j)) for j in range(d)]]
    plt.plot(norms[i+1],label="tour %d"%i, color=colours[i+1])

legend = plt.legend(loc='upper center')
plt.show()
