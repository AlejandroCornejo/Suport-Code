

import sympy
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings

# Output mode to a c++ file
mode = "c"
nnodes = 8
dim = 3
voigt_size = 6
is_symmetric = True # D tangent

w = sympy.Symbol("IntegrationWeight")

B = DefineMatrix('rB', voigt_size, nnodes*dim)
D = DefineMatrix('rD', voigt_size, voigt_size)

if is_symmetric:
    D = DefineSymmetricMatrix('rD', voigt_size, voigt_size)

D = D*w
# A = sympy.simplify(B.transpose()*D)
# C = sympy.simplify(A*B)
A = (B.transpose()*D)
C = (A*B)

out = OutputMatrix_CollectingFactors(C, "rLeftHandSideMatrix", mode, assignment_op="+=")
print(out)
    
# if not is_symmetric:
    # out = OutputMatrix_CollectingFactors(C, "rLeftHandSideMatrix", mode)
    # print(out)
# else:
    # C_sym = sympy.zeros(nnodes*dim,nnodes*dim)
    # for i in range(0, nnodes*dim):
        # for j in range(0, nnodes*dim):
            # if i >= j:
                # C_sym[i,j] = C[i,j]
                # out = OutputScalar_CollectingFactors(C_sym[i,j], "rLeftHandSideMatrix(" + str(i) + "," + str(j) + ")", mode)
                # print(out)
            # else:
                # C_sym[i,j] = 0.0
    # out = OutputMatrix_CollectingFactors(C_sym, "rLeftHandSideMatrix", mode)
    # print(out)