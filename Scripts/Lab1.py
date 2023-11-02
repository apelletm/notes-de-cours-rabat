from fpylll import *
import copy

from fpylll.algorithms.bkz import BKZReduction
from math import log
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
import matplotlib.pyplot as plt

from g6k import *
from g6k.utils.stats import SieveTreeTracer
from g6k.algorithms.bkz import pump_n_jump_bkz_tour

# ReadMe https://readthedocs.org/projects/fpylll/downloads/pdf/latest/
"""
LLL
"""

"""
A = IntegerMatrix.random(30, "qary", k=15, bits=30)
dim = A.nrows


M = GSO.Mat(A, float_type="d") #float_type \in {'d', 'dd', 'qd','mpfr'}
print(M.get_r(0,0))
print(M.d) # lattice dimension
print(M.B) # lattice basis
M.update_gso() # all interesting matrices (GS, mu, ...)
print(M.get_r(1,1))
print(A[0].norm()**2) # squared norm of the first vector
print('before:', [M.get_r(i,i) for i in range(dim)]) #all r_ii

L = LLL.Reduction(M, delta=0.99, eta=0.501) # create new LLL object
L() # run LLL __call__
print('after:', [M.get_r(i,i) for i in range(dim)])
print(A)


A = IntegerMatrix.random(100, "qary", k=50, bits=30)
dim = A.nrows
Ac = copy.deepcopy(A) # Copy A
U = IntegerMatrix.identity(dim)
UinvT = IntegerMatrix.identity(dim)

M_transf = GSO.Mat(Ac, float_type="d", U = U, UinvT = UinvT) #U --  transformation matrix from A to LLL reduced Ac
print(M_transf.inverse_transform_enabled) # LLL вычисляет в процессе U

L = LLL.Reduction(M_transf, delta=0.99, eta=0.501)
L()
print(A[1])
print(Ac[1])
print((U*A)[1])
"""


"""
BKZ
"""

"""
A = IntegerMatrix.random(60, "qary", k=30, bits=30)

bkz = BKZReduction(A) #instance of the class BKZReduction
print('before BKZ:', log(A[0].norm(),2))
bkz(BKZ.EasyParam(20, max_loops=8)) #blocksize 15; maxloops = number of bkz tours
print('after BKZ:', log(A[0].norm(),2))
bkz(BKZ.EasyParam(15, max_loops=8), tracer=True)
print(bkz.trace.get(("tour", 1)).report()) #for tracer=True



k = 40
flags = BKZ.AUTO_ABORT
#print('flags:', flags, BKZ.AUTO_ABORT, BKZ.MAX_LOOPS,BKZ.VERBOSE)
par = BKZ.Param(k, strategies=BKZ.DEFAULT_STRATEGY, max_loops=8, flags=flags) # инстанция класса с пар-ми BKZ
#bkz = BKZ2(A) #
#bkz = BKZ2(GSO.Mat(A)) #
bkz = BKZ2(LLL.Reduction(GSO.Mat(A)))
_ = bkz(par) # run BKZ reduction on A with parameters par
print(log(A[0].norm(),2))


A = IntegerMatrix.random(80, "qary", k=40, bits=45)

LLL.reduction(A)
M = GSO.Mat(A)
M.update_gso()
d = A.nrows

colours = ["#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", "#F15854"]
norms = [[log(M.get_r(i,i)) for i in range(d)]]
plt.plot(norms[0],label="lll", color=colours[0])

beta = 42
tours = 8
par = BKZ.Param(block_size=beta,strategies=BKZ.DEFAULT_STRATEGY)
bkz = BKZ2(M)

for i in range(tours):
    bkz.tour(par)
    norms += [[log(M.get_r(j,j)) for j in range(d)]]
    plt.plot(norms[i+1],label="tour %d"%i, color=colours[i+1])

legend = plt.legend(loc='upper center')
plt.show()
"""


"""
Enumeration
"""
"""
A = IntegerMatrix.random(45, "qary", k=25, bits=30)
M = GSO.Mat(A)
L = LLL.Reduction(M, delta=0.99, eta=0.501)
L()
M.update_gso()

enum = Enumeration(M, strategy=EvaluatorStrategy.BEST_N_SOLUTIONS, sub_solutions=True) #instance of the class Enumeration
res = enum.enumerate(0, 45, 0.9*M.get_r(0, 0), 0) # run enumeration
print(([int(res[0][1][i]) for i in range(len(res[0][1]))]))
for a,b in enum.sub_solutions:
    print(a, b)

print(enum.get_nodes()) # number of visited nodes
"""
