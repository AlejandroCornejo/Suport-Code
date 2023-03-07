

import sympy
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings

# Output mode to a c++ file
mode = "c"
nnodes = 8
dim = 3
voigt_size = 6

B = DefineMatrix('B', voigt_size, nnodes*dim)
D = DefineSymmetricMatrix('D', voigt_size, voigt_size)
A = sympy.simplify(B.transpose()*D)
C = sympy.simplify(A*B)
out = OutputMatrix_CollectingFactors(C, "C", mode)

print(out)