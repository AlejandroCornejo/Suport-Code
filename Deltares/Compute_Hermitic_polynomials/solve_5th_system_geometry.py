import sympy as sp

'''
In this script I am solving the required system of equations to obtain the alpha parameters
such that:
    Pol(x) = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5
    
    dPol_dx (x) = a1 + 2*a2*x + 3*a3*x**2 + 4*a4*x**3 + 5*a5*x**4

fulfiling the boundary conditions at the 3 nodes (position and angle)


Author: Alejandro Cornejo, CIMNE, 2024
'''

x1, y1, alpha1, x2, y2, alpha2, x3, y3, alpha3, a0, a1, a2, a3, a4, a5 = sp.symbols('x1_, y1_, alpha1, x2_, y2_, alpha2, x3_, y3_, alpha3, a0, a1, a2, a3, a4, a5')

# position BC
eq1 = sp.Eq(a0 + a1*x1 + a2*x1**2 + a3*x1**3 + a4*x1**4 + a5*x1**5, y1)
eq2 = sp.Eq(a0 + a1*x2 + a2*x2**2 + a3*x2**3 + a4*x2**4 + a5*x2**5, y2)
eq3 = sp.Eq(a0 + a1*x3 + a2*x3**2 + a3*x3**3 + a4*x3**4 + a5*x3**5, y3)

# angle BC
eq4 = sp.Eq(a1 + 2*a2*x1 + 3*a3*x1**2 + 4*a4*x1**3 + 5*a5*x1**4, alpha1)
eq5 = sp.Eq(a1 + 2*a2*x2 + 3*a3*x2**2 + 4*a4*x2**3 + 5*a5*x2**4, alpha2)
eq6 = sp.Eq(a1 + 2*a2*x3 + 3*a3*x3**2 + 4*a4*x3**3 + 5*a5*x3**4, alpha3)

solution = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6), (a0, a1, a2, a3, a4, a5))

#print("The value of alpha_0 is: ", solution[a0])
#print("The value of alpha_1 is: ", solution[a1])
#print("The value of alpha_2 is: ", solution[a2])
#print("The value of alpha_3 is: ", solution[a3])
#print("The value of alpha_4 is: ", solution[a4])
print("The value of alpha_5 is: ", solution[a5])