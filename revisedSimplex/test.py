###test.py for optimization#########
from se import *
from rsimplex import *

a = SimplexElement()
a.readFromFile("ques.txt")
solver = RevisedSimplex(a)
solver.solve()
print("Optimized Value  =", solver.optimizedVal)
