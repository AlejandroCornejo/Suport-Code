import sympy as sp

'''
In this script I am solving the required system of equations to obtain the alpha parameters
such that:
    v = alpha_0 + alpha_1*xi + alpha_2*xi**2 + alpha_3*xi**3

fulfiling the boundary conditions at the 2 nodes

Euler-Bernoulli theory: theta = dv/dx
'''

v1, theta1, v2, theta2, L, alpha0, alpha1, alpha2, alpha3 = sp.symbols('v1 theta1 v2 theta2 L alpha0 alpha1 alpha2 alpha3')

eq1 = sp.Eq(alpha0 - alpha1 + alpha2 - alpha3, v1)
eq2 = sp.Eq(alpha0 + alpha1 + alpha2 + alpha3, v2)

eq3 = sp.Eq(2.0 / L * (alpha1 - 2.0*alpha2 + 3.0 * alpha3), theta1)
eq4 = sp.Eq(2.0 / L * (alpha1 + 2.0*alpha2 + 3.0 * alpha3), theta2)

solution = sp.solve((eq1, eq2, eq3, eq4), (alpha0, alpha1, alpha2, alpha3))

print("The value of alpha_0 is: ", solution[alpha0])
print("The value of alpha_1 is: ", solution[alpha1])
print("The value of alpha_2 is: ", solution[alpha2])
print("The value of alpha_3 is: ", solution[alpha3])

'''
Output:
    The value of alpha_0 is:  0.125*L*theta1 - 0.125*L*theta2 + 0.5*v1 + 0.5*v2
    The value of alpha_1 is:  -0.125*L*theta1 - 0.125*L*theta2 - 0.75*v1 + 0.75*v2
    The value of alpha_2 is:  -0.125*L*theta1 + 0.125*L*theta2
    The value of alpha_3 is:  0.125*L*theta1 + 0.125*L*theta2 + 0.25*v1 - 0.25*v2 
'''

# let's check...
L = 1.0
v1 = 0.1
v2 = 0.25
theta1 = -0.0025
theta2 = 0.001

alpha_0 =  0.125*L*theta1 - 0.125*L*theta2 + 0.5*v1 + 0.5*v2
alpha_1 =  -0.125*L*theta1 - 0.125*L*theta2 - 0.75*v1 + 0.75*v2
alpha_2 =  -0.125*L*theta1 + 0.125*L*theta2
alpha_3 =  0.125*L*theta1 + 0.125*L*theta2 + 0.25*v1 - 0.25*v2 

def ComputeV(alpha_0, alpha_1, alpha_2, alpha_3, xi):
    return alpha_0 + alpha_1*xi + alpha_2*xi**2 + alpha_3*xi**3
def ComputeTheta(alpha_0, alpha_1, alpha_2, alpha_3, xi):
    return alpha_1 + 2*alpha_2*xi + 3*alpha_3*xi**2

print("The calculated v at xi=-1 is: ", ComputeV(alpha_0, alpha_1, alpha_2, alpha_3, -1))
print("The calculated v at xi= 1 is: ", ComputeV(alpha_0, alpha_1, alpha_2, alpha_3,  1))
print("The calculated Theta at xi=-1 is: ", 2/L*ComputeTheta(alpha_0, alpha_1, alpha_2, alpha_3, -1))