
import sympy as sp
import numpy as np

'''
In this script I am solving the required system of equations to obtain the "a" parameters
of the deflection polynomial of a curved beam

The curved beam is parametrized by a previously known 5th order polynomial:

    y(x)   = p0 + p1*x + p2*x**2 + p3*x**3 + p4*x**4 + p5*x**5               (1)
    dy/dx  = p1 + 2*p2*x + 3*p3*x**2 + 4*p4*x**3 + 5*p5*x**4                 (2)
    k0(x)  = 2*p2 + 6*p3*x + 12*p4*x**2 + 20*p5*x**3                         (3)

    x = 0.5*(x2-x1)XI + 0.5*(x1+x2)

such that:
                              |--> xi
Isoparametric element:   O----O----O
                        1    3     2

    v = a0 + a1*xi + a2*xi**2 + a3*xi**3 + a4*xi**4 + a5*xi**5               (4)

    dv_dxi   = a1 + 2*a2*xi + 3*a3*xi**2 + 4*a4*xi**3 + 5*a5*xi**4
    d2v_dxi2 = 2*a2 + 6*a3*xi + 12*a4*xi**2 + 20*a5*xi**3
    d3v_dxi3 = 6*a3 + 24*a4*xi + 60*a5*xi**2

fulfiling the boundary conditions at the 3 nodes:

    - deflection at each node, i.e.             v(xi=-1) = v1 ; v(xi=0) = v3  ;  v(xi=1) = v2      Eqs (5,6,7)

    - total rotation fulfilled at each node:    theta(xi=-1) = theta1 ; and so on.                 Eqs (8,9,10)

Curved Timoshenko theory:    theta = dv/ds   +   (phi * L**2 / 12.0) * d3v/ds3   +   k0 * u 

Reminder:    dN/ds = (1 / J) * dN/dxi

being        J = 0.5 * (x2-x1) * sqrt(1+(dy/dx)**2)

phi is the shear parameter

Author: Alejandro Cornejo, CIMNE, 2024

'''

# Let's define all symbols:

# Coordinates of the nodes, known
# x1, x2, x3, y1, y2, y3 = sp.symbols('x1 x2 x3 y1 y2 y3')

# isoparametric coordinate xi
xi = sp.symbols('xi')

# components of the known parametrization of the geometry
# p0, p1, p2, p3, p4, p5 = sp.symbols('p0 p1 p2 p3 p4 p5')

# shear parameter
k_s = sp.symbols('k_s') # k_s = E·I / G·As

'''
# Let's compute the dy_dx at each node
dydx1 = p1 + 2.0 * p2 *x1 + 3.0 * p3 * x1**2 + 4.0 * p4 * x1**3 + 5.0 * p5 * x1**4
dydx2 = p1 + 2.0 * p2 *x2 + 3.0 * p3 * x2**2 + 4.0 * p4 * x2**3 + 5.0 * p5 * x2**4
dydx3 = p1 + 2.0 * p2 *x3 + 3.0 * p3 * x3**2 + 4.0 * p4 * x3**3 + 5.0 * p5 * x3**4

# IDEM with initial curvatures k0
k01 = 2.0 * p2 + 6.0 * p3 * x1 + 12.0 * p4 * x1**2 + 20.0 * p5 * x1**3
k02 = 2.0 * p2 + 6.0 * p3 * x2 + 12.0 * p4 * x2**2 + 20.0 * p5 * x2**3
k03 = 2.0 * p2 + 6.0 * p3 * x3 + 12.0 * p4 * x3**2 + 20.0 * p5 * x3**3

# We computhe the Jacobian at each node
J1 = 0.5 * (x2-x1) * sp.sqrt(1.0 + dydx1**2)
J2 = 0.5 * (x2-x1) * sp.sqrt(1.0 + dydx2**2)
J3 = 0.5 * (x2-x1) * sp.sqrt(1.0 + dydx3**2)
'''

k01, k02, k03, J1, J2, J3 = sp.symbols('k01 k02 k03 J1 J2 J3')

# The nodal unknowns
u1, v1, theta1, u2, v2, theta2, u3, v3, theta3 = sp.symbols('u1 v1 theta1 u2 v2 theta2 u3 v3 theta3')

# The components to be calculated of the deflection polynomial
a0, a1, a2, a3, a4, a5 = sp.symbols('a0 a1 a2 a3 a4 a5')

# Enforcing deflections...
eq1 = sp.Eq(a0 - a1 + a2 - a3 + a4 - a5, v1)
eq2 = sp.Eq(a0 + a1 + a2 + a3 + a4 + a5, v2)
eq3 = sp.Eq(a0                         , v3)

# Enforcing total rotations...
eq4 = sp.Eq((1.0 / J1) * (a1 - 2.0 * a2 + 3.0 * a3 - 4.0 * a4 + 5.0 * a5) + k_s * (1.0 / J1**3) * (6.0 * a3 - 24.0 * a4 + 60.0 * a5) + u1 * k01, theta1)
eq5 = sp.Eq((1.0 / J2) * (a1 + 2.0 * a2 + 3.0 * a3 + 4.0 * a4 + 5.0 * a5) + k_s * (1.0 / J2**3) * (6.0 * a3 + 24.0 * a4 + 60.0 * a5) + u2 * k02, theta2)
eq6 = sp.Eq((1.0 / J3) * (a1)                                             + k_s * (1.0 / J3**3) * (6.0 * a3)                         + u3 * k03, theta3)

solution = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6), (a0, a1, a2, a3, a4, a5))

# print("The value of a0 is: ", solution[a0])
# print("The value of a1 is: ", solution[a1])
# print("The value of a2 is: ", solution[a2])
# print("The value of a3 is: ", solution[a3])
# print("The value of a4 is: ", solution[a4])
# print("The value of a5 is: ", solution[a5])



'''
Output:

--
The value of a0 is:  

v3

--
The value of a1 is:  (multiplied by xi)

(-3.0*J1**3*J2**2*k01*k_s*u1 + 3.0*J1**3*J2**2*k_s*theta1 - 36.0*J1**3*k01*k_s**2*u1 + 36.0*J1**3*k_s**2*theta1 - 3.0*J1**2*J2**3*k02*k_s*u2 + 3.0*J1**2*J2**3*k_s*theta2 - 2.0*J1**2*J2**2*J3**3*k03*u3 + 2.0*J1**2*J2**2*J3**3*theta3 + 15.0*J1**2*J2**2*k_s*v1 - 15.0*J1**2*J2**2*k_s*v2 - 39.0*J1**2*J3**3*k03*k_s*u3 + 39.0*J1**2*J3**3*k_s*theta3 + 216.0*J1**2*k_s**2*v1 - 144.0*J1**2*k_s**2*v2 - 72.0*J1**2*k_s**2*v3 - 36.0*J2**3*k02*k_s**2*u2 + 36.0*J2**3*k_s**2*theta2 - 39.0*J2**2*J3**3*k03*k_s*u3 + 39.0*J2**2*J3**3*k_s*theta3 + 144.0*J2**2*k_s**2*v1 - 216.0*J2**2*k_s**2*v2 + 72.0*J2**2*k_s**2*v3 - 648.0*J3**3*k03*k_s**2*u3 + 648.0*J3**3*k_s**2*theta3 + 2160.0*k_s**3*v1 - 2160.0*k_s**3*v2)/(2.0*J1**2*J2**2*J3**2 - 24.0*J1**2*J2**2*k_s + 39.0*J1**2*J3**2*k_s - 324.0*J1**2*k_s**2 + 39.0*J2**2*J3**2*k_s - 324.0*J2**2*k_s**2 + 648.0*J3**2*k_s**2 - 4320.0*k_s**3)

--
The value of a2 is:  (multiplied by xi**2)

(-2.0*J1**3*J2**2*J3**2*k01*u1 + 2.0*J1**3*J2**2*J3**2*theta1 + 24.0*J1**3*J2**2*k01*k_s*u1 - 24.0*J1**3*J2**2*k_s*theta1 - 54.0*J1**3*J3**2*k01*k_s*u1 + 54.0*J1**3*J3**2*k_s*theta1 + 
360.0*J1**3*k01*k_s**2*u1 - 360.0*J1**3*k_s**2*theta1 + 2.0*J1**2*J2**3*J3**2*k02*u2 - 2.0*J1**2*J2**3*J3**2*theta2 - 24.0*J1**2*J2**3*k02*k_s*u2 + 24.0*J1**2*J2**3*k_s*theta2 + 8.0*J1**2*J2**2*J3**2*v1 + 
8.0*J1**2*J2**2*J3**2*v2 - 16.0*J1**2*J2**2*J3**2*v3 - 96.0*J1**2*J2**2*k_s*v1 - 96.0*J1**2*J2**2*k_s*v2 + 192.0*J1**2*J2**2*k_s*v3 - 96.0*J1**2*J3**3*k03*k_s*u3 + 96.0*J1**2*J3**3*k_s*theta3 + 207.0*J1**2*J3**2*k_s*v1 + 57.0*J1**2*J3**2*k_s*v2 - 264.0*J1**2*J3**2*k_s*v3 - 1188.0*J1**2*k_s**2*v1 - 828.0*J1**2*k_s**2*v2 + 2016.0*J1**2*k_s**2*v3 + 54.0*J2**3*J3**2*k02*k_s*u2 - 54.0*J2**3*J3**2*k_s*theta2 - 360.0*J2**3*k02*k_s**2*u2 + 360.0*J2**3*k_s**2*theta2 + 96.0*J2**2*J3**3*k03*k_s*u3 - 96.0*J2**2*J3**3*k_s*theta3 + 57.0*J2**2*J3**2*k_s*v1 + 207.0*J2**2*J3**2*k_s*v2 - 264.0*J2**2*J3**2*k_s*v3 - 828.0*J2**2*k_s**2*v1 - 1188.0*J2**2*k_s**2*v2 + 2016.0*J2**2*k_s**2*v3 + 1296.0*J3**2*k_s**2*v1 + 1296.0*J3**2*k_s**2*v2 - 2592.0*J3**2*k_s**2*v3 - 8640.0*k_s**3*v1 - 8640.0*k_s**3*v2 + 17280.0*k_s**3*v3)/(8.0*J1**2*J2**2*J3**2 - 96.0*J1**2*J2**2*k_s + 156.0*J1**2*J3**2*k_s - 1296.0*J1**2*k_s**2 + 156.0*J2**2*J3**2*k_s - 1296.0*J2**2*k_s**2 + 2592.0*J3**2*k_s**2 - 17280.0*k_s**3)

--
The value of a3 is:  (multiplied by xi**3)

(J1**3*J2**2*J3**2*k01*u1 - J1**3*J2**2*J3**2*theta1 + 12.0*J1**3*J3**2*k01*k_s*u1 - 12.0*J1**3*J3**2*k_s*theta1 + J1**2*J2**3*J3**2*k02*u2 - J1**2*J2**3*J3**2*theta2 + 8.0*J1**2*J2**2*J3**3*k03*u3 - 8.0*J1**2*J2**2*J3**3*theta3 - 5.0*J1**2*J2**2*J3**2*v1 + 5.0*J1**2*J2**2*J3**2*v2 + 108.0*J1**2*J3**3*k03*k_s*u3 - 108.0*J1**2*J3**3*k_s*theta3 - 72.0*J1**2*J3**2*k_s*v1 + 48.0*J1**2*J3**2*k_s*v2 + 24.0*J1**2*J3**2*k_s*v3 + 12.0*J2**3*J3**2*k02*k_s*u2 - 12.0*J2**3*J3**2*k_s*theta2 + 108.0*J2**2*J3**3*k03*k_s*u3 - 108.0*J2**2*J3**3*k_s*theta3 - 48.0*J2**2*J3**2*k_s*v1 + 72.0*J2**2*J3**2*k_s*v2 - 24.0*J2**2*J3**2*k_s*v3 + 1440.0*J3**3*k03*k_s**2*u3 - 1440.0*J3**3*k_s**2*theta3 - 720.0*J3**2*k_s**2*v1 + 720.0*J3**2*k_s**2*v2)/(4.0*J1**2*J2**2*J3**2 - 48.0*J1**2*J2**2*k_s + 78.0*J1**2*J3**2*k_s 
- 648.0*J1**2*k_s**2 + 78.0*J2**2*J3**2*k_s - 648.0*J2**2*k_s**2 + 1296.0*J3**2*k_s**2 - 8640.0*k_s**3)

--
The value of a4 is:  (multiplied by xi**4)

(2.0*J1**3*J2**2*J3**2*k01*u1 - 2.0*J1**3*J2**2*J3**2*theta1 - 24.0*J1**3*J2**2*k01*k_s*u1 + 24.0*J1**3*J2**2*k_s*theta1 + 54.0*J1**3*J3**2*k01*k_s*u1 - 54.0*J1**3*J3**2*k_s*theta1 - 360.0*J1**3*k01*k_s**2*u1 + 360.0*J1**3*k_s**2*theta1 - 2.0*J1**2*J2**3*J3**2*k02*u2 + 2.0*J1**2*J2**3*J3**2*theta2 + 24.0*J1**2*J2**3*k02*k_s*u2 - 24.0*J1**2*J2**3*k_s*theta2 - 4.0*J1**2*J2**2*J3**2*v1 - 4.0*J1**2*J2**2*J3**2*v2 + 8.0*J1**2*J2**2*J3**2*v3 + 48.0*J1**2*J2**2*k_s*v1 + 48.0*J1**2*J2**2*k_s*v2 - 96.0*J1**2*J2**2*k_s*v3 + 96.0*J1**2*J3**3*k03*k_s*u3 - 96.0*J1**2*J3**3*k_s*theta3 - 129.0*J1**2*J3**2*k_s*v1 + 21.0*J1**2*J3**2*k_s*v2 + 108.0*J1**2*J3**2*k_s*v3 + 540.0*J1**2*k_s**2*v1 + 180.0*J1**2*k_s**2*v2 - 720.0*J1**2*k_s**2*v3 - 54.0*J2**3*J3**2*k02*k_s*u2 + 54.0*J2**3*J3**2*k_s*theta2 + 360.0*J2**3*k02*k_s**2*u2 - 360.0*J2**3*k_s**2*theta2 - 96.0*J2**2*J3**3*k03*k_s*u3 + 96.0*J2**2*J3**3*k_s*theta3 + 21.0*J2**2*J3**2*k_s*v1 - 129.0*J2**2*J3**2*k_s*v2 + 108.0*J2**2*J3**2*k_s*v3 + 180.0*J2**2*k_s**2*v1 + 540.0*J2**2*k_s**2*v2 - 720.0*J2**2*k_s**2*v3)/(8.0*J1**2*J2**2*J3**2 - 96.0*J1**2*J2**2*k_s + 156.0*J1**2*J3**2*k_s - 1296.0*J1**2*k_s**2 + 156.0*J2**2*J3**2*k_s - 1296.0*J2**2*k_s**2 + 2592.0*J3**2*k_s**2 - 17280.0*k_s**3)

--
The value of a5 is:  (multiplied by xi**5)

(-J1**3*J2**2*J3**2*k01*u1 + J1**3*J2**2*J3**2*theta1 + 6.0*J1**3*J2**2*k01*k_s*u1 - 6.0*J1**3*J2**2*k_s*theta1 - 12.0*J1**3*J3**2*k01*k_s*u1 + 12.0*J1**3*J3**2*k_s*theta1 + 72.0*J1**3*k01*k_s**2*u1 - 72.0*J1**3*k_s**2*theta1 - J1**2*J2**3*J3**2*k02*u2 + J1**2*J2**3*J3**2*theta2 + 6.0*J1**2*J2**3*k02*k_s*u2 - 6.0*J1**2*J2**3*k_s*theta2 - 4.0*J1**2*J2**2*J3**3*k03*u3 + 4.0*J1**2*J2**2*J3**3*theta3 + 3.0*J1**2*J2**2*J3**2*v1 - 3.0*J1**2*J2**2*J3**2*v2 - 6.0*J1**2*J2**2*k_s*v1 + 6.0*J1**2*J2**2*k_s*v2 - 30.0*J1**2*J3**3*k03*k_s*u3 + 30.0*J1**2*J3**3*k_s*theta3 + 33.0*J1**2*J3**2*k_s*v1 - 9.0*J1**2*J3**2*k_s*v2 - 24.0*J1**2*J3**2*k_s*v3 - 108.0*J1**2*k_s**2*v1 - 36.0*J1**2*k_s**2*v2 + 144.0*J1**2*k_s**2*v3 - 12.0*J2**3*J3**2*k02*k_s*u2 + 12.0*J2**3*J3**2*k_s*theta2 + 72.0*J2**3*k02*k_s**2*u2 
- 72.0*J2**3*k_s**2*theta2 - 30.0*J2**2*J3**3*k03*k_s*u3 + 30.0*J2**2*J3**3*k_s*theta3 + 9.0*J2**2*J3**2*k_s*v1 - 33.0*J2**2*J3**2*k_s*v2 + 24.0*J2**2*J3**2*k_s*v3 + 36.0*J2**2*k_s**2*v1 + 108.0*J2**2*k_s**2*v2 - 144.0*J2**2*k_s**2*v3 - 144.0*J3**3*k03*k_s**2*u3 + 144.0*J3**3*k_s**2*theta3 + 72.0*J3**2*k_s**2*v1 - 72.0*J3**2*k_s**2*v2)/(4.0*J1**2*J2**2*J3**2 - 48.0*J1**2*J2**2*k_s + 78.0*J1**2*J3**2*k_s - 648.0*J1**2*k_s**2 + 78.0*J2**2*J3**2*k_s - 648.0*J2**2*k_s**2 + 1296.0*J3**2*k_s**2 - 8640.0*k_s**3)

'''


'''
Now we want to transform the:

    v = a0 + a1*xi + a2*xi**2 + a3*xi**3 + a4*xi**4 + a5*xi**5

to be

    v = {Nu1, Nv1, Nt1, Nu2, Nv2, Nt2, Nu3, Nv3, Nt3} · {u1, v1, theta1, u2, v2, theta2, u3, v3, theta3} = N·u

------------------

'''

Nu1, Nv1, Nt1, Nu2, Nv2, Nt2, Nu3, Nv3, Nt3 = sp.symbols('Nu1 Nv1 Nt1 Nu2 Nv2 Nt2 Nu3 Nv3 Nt3')

poly = solution[a0] + solution[a1]*xi +solution[a2]*xi**2 + solution[a3]*xi**3 + solution[a4]*xi**4 + solution[a5]*xi**5

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 1.0)
# eq_n1  = sp.Eq(Nu1, poly)
# solution_Nu1 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nu1, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nu1: \n", sp.ccode(solution_Nu1[Nu1]))
# print("dNu1_dxi: \n", sp.diff(solution_Nu1[Nu1], xi))
# print("d2Nu1_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nu1[Nu1], xi), xi)))
# print("d3Nu1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nu1[Nu1], xi), xi), xi)))
# print("d4Nu1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nu1[Nu1], xi), xi), xi), xi)))

'''
Nu1:

Nu1:
 (-pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*pow(xi, 5) + pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*pow(xi, 4) + pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*pow(xi, 3) - pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*pow(xi, 2) + 6.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*pow(xi, 5) - 12.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*pow(xi, 4) + 12.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*pow(xi, 2) - 6.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*xi - 12.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*pow(xi, 5) + 27.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*pow(xi, 4) + 12.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*pow(xi, 3) - 27.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*pow(xi, 2) + 72.0*pow(J1, 3)*k01*pow(k_s, 2)*pow(xi, 5) - 180.0*pow(J1, 3)*k01*pow(k_s, 2)*pow(xi, 4) + 180.0*pow(J1, 3)*k01*pow(k_s, 2)*pow(xi, 2) - 72.0*pow(J1, 3)*k01*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

dNu1_dxi: 
 (-5*J1_**3*J2_**2*J3_**2*k01_*xi**4 + 4*J1_**3*J2_**2*J3_**2*k01_*xi**3 + 3*J1_**3*J2_**2*J3_**2*k01_*xi**2 - 2*J1_**3*J2_**2*J3_**2*k01_*xi + 30.0*J1_**3*J2_**2*k01_*k_s*xi**4 
- 48.0*J1_**3*J2_**2*k01_*k_s*xi**3 + 24.0*J1_**3*J2_**2*k01_*k_s*xi - 6.0*J1_**3*J2_**2*k01_*k_s - 60.0*J1_**3*J3_**2*k01_*k_s*xi**4 + 108.0*J1_**3*J3_**2*k01_*k_s*xi**3 + 36.0*J1_**3*J3_**2*k01_*k_s*xi**2 - 54.0*J1_**3*J3_**2*k01_*k_s*xi + 360.0*J1_**3*k01_*k_s**2*xi**4 - 720.0*J1_**3*k01_*k_s**2*xi**3 + 360.0*J1_**3*k01_*k_s**2*xi - 72.0*J1_**3*k01_*k_s**2)/(4.0*J1_**2*J2_**2*J3_**2 - 48.0*J1_**2*J2_**2*k_s + 78.0*J1_**2*J3_**2*k_s - 648.0*J1_**2*k_s**2 + 78.0*J2_**2*J3_**2*k_s - 648.0*J2_**2*k_s**2 + 1296.0*J3_**2*k_s**2 - 8640.0*k_s**3)

d2Nu1_dxi: 
 (-20*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*k01_*pow(xi, 3) + 12*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*k01_*pow(xi, 2) + 6*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*k01_*xi - 2*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*k01_ + 120.0*pow(J1_, 3)*pow(J2_, 2)*k01_*k_s*pow(xi, 3) - 144.0*pow(J1_, 3)*pow(J2_, 2)*k01_*k_s*pow(xi, 2) + 24.0*pow(J1_, 3)*pow(J2_, 2)*k01_*k_s - 240.0*pow(J1_, 3)*pow(J3_, 2)*k01_*k_s*pow(xi, 3) + 324.0*pow(J1_, 3)*pow(J3_, 2)*k01_*k_s*pow(xi, 2) + 72.0*pow(J1_, 3)*pow(J3_, 2)*k01_*k_s*xi - 54.0*pow(J1_, 3)*pow(J3_, 2)*k01_*k_s + 1440.0*pow(J1_, 3)*k01_*pow(k_s, 2)*pow(xi, 3) - 2160.0*pow(J1_, 3)*k01_*pow(k_s, 2)*pow(xi, 2) + 360.0*pow(J1_, 3)*k01_*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d3Nu1_dxi: 
 (-60*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*pow(xi, 2) + 24*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*xi + 6*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01 + 360.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*pow(xi, 2) - 288.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*xi - 720.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*pow(xi, 2) + 648.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*xi + 72.0*pow(J1, 3)*pow(J3, 2)*k01*k_s + 4320.0*pow(J1, 3)*k01*pow(k_s, 2)*pow(xi, 2) - 4320.0*pow(J1, 3)*k01*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))


d4Nu1_dxi:
 (-120*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01*xi + 24*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*k01 + 720.0*pow(J1, 3)*pow(J2, 2)*k01*k_s*xi - 288.0*pow(J1, 3)*pow(J2, 2)*k01*k_s - 1440.0*pow(J1, 3)*pow(J3, 2)*k01*k_s*xi + 648.0*pow(J1, 3)*pow(J3, 2)*k01*k_s + 8640.0*pow(J1, 3)*k01*pow(k_s, 2)*xi - 4320.0*pow(J1, 3)*k01*pow(k_s, 2))/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))
'''

# eq_n2  = sp.Eq(    v1, 1.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nv1, poly)
# solution_Nv1 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nv1, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nv1: \n", sp.ccode(solution_Nv1[Nv1]))
# print("dNv1_dxi: \n", sp.ccode(sp.diff(solution_Nv1[Nv1], xi)))
# print("d2Nv1_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nv1[Nv1], xi), xi)))
# print("d3Nv1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nv1[Nv1], xi), xi), xi)))
# print("\n d4Nv1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nv1[Nv1], xi), xi), xi), xi)))

'''
Nv1:

(6.0*J1_**2*J2_**2*J3_**2*xi**5 - 4.0*J1_**2*J2_**2*J3_**2*xi**4 - 10.0*J1_**2*J2_**2*J3_**2*xi**3 + 8.0*J1_**2*J2_**2*J3_**2*xi**2 - 12.0*J1_**2*J2_**2*k_s*xi**5 + 48.0*J1_**2*J2_**2*k_s*xi**4 - 96.0*J1_**2*J2_**2*k_s*xi**2 + 60.0*J1_**2*J2_**2*k_s*xi + 66.0*J1_**2*J3_**2*k_s*xi**5 - 129.0*J1_**2*J3_**2*k_s*xi**4 - 144.0*J1_**2*J3_**2*k_s*xi**3 + 207.0*J1_**2*J3_**2*k_s*xi**2 - 216.0*J1_**2*k_s**2*xi**5 + 540.0*J1_**2*k_s**2*xi**4 - 1188.0*J1_**2*k_s**2*xi**2 + 864.0*J1_**2*k_s**2*xi + 18.0*J2_**2*J3_**2*k_s*xi**5 + 21.0*J2_**2*J3_**2*k_s*xi**4 - 96.0*J2_**2*J3_**2*k_s*xi**3 + 57.0*J2_**2*J3_**2*k_s*xi**2 
+ 72.0*J2_**2*k_s**2*xi**5 + 180.0*J2_**2*k_s**2*xi**4 - 828.0*J2_**2*k_s**2*xi**2 + 576.0*J2_**2*k_s**2*xi + 144.0*J3_**2*k_s**2*xi**5 - 1440.0*J3_**2*k_s**2*xi**3 + 1296.0*J3_**2*k_s**2*xi**2 - 8640.0*k_s**3*xi**2 + 8640.0*k_s**3*xi)/(8.0*J1_**2*J2_**2*J3_**2 - 96.0*J1_**2*J2_**2*k_s + 156.0*J1_**2*J3_**2*k_s - 1296.0*J1_**2*k_s**2 + 156.0*J2_**2*J3_**2*k_s - 1296.0*J2_**2*k_s**2 + 2592.0*J3_**2*k_s**2 - 
17280.0*k_s**3)

dNv1_dxi: 
 (30.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 4) - 16.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 30.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) + 16.0*pow(J1_, 
2)*pow(J2_, 2)*pow(J3_, 2)*xi - 60.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 4) + 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 3) - 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s*xi + 60.0*pow(J1_, 
2)*pow(J2_, 2)*k_s + 330.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) - 516.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) - 432.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 414.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi - 1080.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 4) + 2160.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) - 2376.0*pow(J1_, 2)*pow(k_s, 2)*xi + 864.0*pow(J1_, 2)*pow(k_s, 2) + 90.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) + 84.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) - 288.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 114.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi + 360.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 4) + 720.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) - 1656.0*pow(J2_, 2)*pow(k_s, 2)*xi + 576.0*pow(J2_, 2)*pow(k_s, 2) + 720.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 4) - 4320.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 2) + 2592.0*pow(J3_, 2)*pow(k_s, 2)*xi - 17280.0*pow(k_s, 3)*xi + 8640.0*pow(k_s, 3))/(8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 156.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J1_, 2)*pow(k_s, 2) + 156.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J2_, 2)*pow(k_s, 
2) + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

d2Nv1_dxi: 
 (120.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 48.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) - 60.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*xi + 16.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 240.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 3) + 576.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 2) - 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 1320.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) - 1548.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) - 864.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi + 414.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 4320.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) + 6480.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 2) - 2376.0*pow(J1_, 2)*pow(k_s, 2) + 360.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 252.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) - 576.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi + 114.0*pow(J2_, 2)*pow(J3_, 2)*k_s + 1440.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) + 2160.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 2) - 1656.0*pow(J2_, 2)*pow(k_s, 2) + 2880.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 3) - 8640.0*pow(J3_, 2)*pow(k_s, 2)*xi + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))/(8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 156.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J1_, 2)*pow(k_s, 2) + 156.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J2_, 2)*pow(k_s, 2) + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

 d3Nv1_dxi: 
 (360.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*pow(xi, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*xi - 60.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 720.0*pow(J1, 2)*pow(J2, 2)*k_s*pow(xi, 2) + 1152.0*pow(J1, 2)*pow(J2, 2)*k_s*xi + 3960.0*pow(J1, 2)*pow(J3, 2)*k_s*pow(xi, 2) - 3096.0*pow(J1, 2)*pow(J3, 2)*k_s*xi - 864.0*pow(J1, 2)*pow(J3, 2)*k_s - 12960.0*pow(J1, 2)*pow(k_s, 2)*pow(xi, 2) + 12960.0*pow(J1, 2)*pow(k_s, 2)*xi + 1080.0*pow(J2, 2)*pow(J3, 2)*k_s*pow(xi, 2) + 504.0*pow(J2, 2)*pow(J3, 2)*k_s*xi - 576.0*pow(J2, 2)*pow(J3, 2)*k_s + 4320.0*pow(J2, 2)*pow(k_s, 2)*pow(xi, 2) + 4320.0*pow(J2, 2)*pow(k_s, 2)*xi + 8640.0*pow(J3, 2)*pow(k_s, 2)*pow(xi, 2) - 8640.0*pow(J3, 2)*pow(k_s, 2))/(8.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*k_s + 156.0*pow(J1, 2)*pow(J3, 2)*k_s - 1296.0*pow(J1, 2)*pow(k_s, 2) + 156.0*pow(J2, 2)*pow(J3, 2)*k_s - 1296.0*pow(J2, 2)*pow(k_s, 2) 
+ 2592.0*pow(J3, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

 d4Nv1_dxi: 
 (720.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*xi - 96.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 1440.0*pow(J1, 2)*pow(J2, 2)*k_s*xi + 1152.0*pow(J1, 2)*pow(J2, 2)*k_s + 7920.0*pow(J1, 2)*pow(J3, 2)*k_s*xi - 3096.0*pow(J1, 2)*pow(J3, 2)*k_s - 25920.0*pow(J1, 2)*pow(k_s, 2)*xi + 12960.0*pow(J1, 2)*pow(k_s, 2) + 2160.0*pow(J2, 2)*pow(J3, 2)*k_s*xi + 504.0*pow(J2, 2)*pow(J3, 2)*k_s + 8640.0*pow(J2, 2)*pow(k_s, 2)*xi + 4320.0*pow(J2, 2)*pow(k_s, 2) + 17280.0*pow(J3, 2)*pow(k_s, 2)*xi)/(8.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*k_s + 156.0*pow(J1, 2)*pow(J3, 2)*k_s - 1296.0*pow(J1, 2)*pow(k_s, 2) + 156.0*pow(J2, 2)*pow(J3, 2)*k_s - 1296.0*pow(J2, 2)*pow(k_s, 2) + 2592.0*pow(J3, 2)*pow(k_s, 2) - 
17280.0*pow(k_s, 3))

'''

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 1.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nt1, poly)
# solution_Nt1 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nt1, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nt1: \n", sp.ccode(solution_Nt1[Nt1]))
# print("dNt1_dxi: \n", sp.ccode(sp.diff(solution_Nt1[Nt1], xi)))
# print("d2Nt1_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nt1[Nt1], xi), xi)))
# print("d3Nt1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nt1[Nt1], xi), xi), xi)))
# print("\n d4Nt1_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nt1[Nt1], xi), xi), xi), xi)))

'''
Nt1:

(pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*pow(xi, 5) - pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*pow(xi, 4) - pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*pow(xi, 3) + pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*pow(xi, 2) - 6.0*pow(J1, 3)*pow(J2, 2)*k_s*pow(xi, 5) + 12.0*pow(J1, 3)*pow(J2, 2)*k_s*pow(xi, 4) - 12.0*pow(J1, 3)*pow(J2, 2)*k_s*pow(xi, 2) + 6.0*pow(J1, 3)*pow(J2, 2)*k_s*xi + 12.0*pow(J1, 3)*pow(J3, 2)*k_s*pow(xi, 5) - 27.0*pow(J1, 3)*pow(J3, 2)*k_s*pow(xi, 4) - 12.0*pow(J1, 3)*pow(J3, 2)*k_s*pow(xi, 3) + 27.0*pow(J1, 3)*pow(J3, 2)*k_s*pow(xi, 2) - 72.0*pow(J1, 3)*pow(k_s, 2)*pow(xi, 5) + 180.0*pow(J1, 3)*pow(k_s, 2)*pow(xi, 4) - 180.0*pow(J1, 3)*pow(k_s, 2)*pow(xi, 2) + 72.0*pow(J1, 3)*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

dNt1_dxi: 
 (5*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 4) - 4*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 3*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) + 2*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*xi - 
30.0*pow(J1_, 3)*pow(J2_, 2)*k_s*pow(xi, 4) + 48.0*pow(J1_, 3)*pow(J2_, 2)*k_s*pow(xi, 3) - 24.0*pow(J1_, 3)*pow(J2_, 2)*k_s*xi + 6.0*pow(J1_, 3)*pow(J2_, 2)*k_s + 60.0*pow(J1_, 3)*pow(J3_, 2)*k_s*pow(xi, 4) - 108.0*pow(J1_, 3)*pow(J3_, 2)*k_s*pow(xi, 3) - 36.0*pow(J1_, 3)*pow(J3_, 2)*k_s*pow(xi, 2) + 54.0*pow(J1_, 3)*pow(J3_, 2)*k_s*xi - 360.0*pow(J1_, 3)*pow(k_s, 2)*pow(xi, 4) + 720.0*pow(J1_, 3)*pow(k_s, 2)*pow(xi, 3) - 360.0*pow(J1_, 3)*pow(k_s, 2)*xi + 72.0*pow(J1_, 3)*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

d2Nt1_dxi: 
 (20*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 12*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) - 6*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2)*xi + 2*pow(J1_, 3)*pow(J2_, 2)*pow(J3_, 2) - 120.0*pow(J1_, 3)*pow(J2_, 2)*k_s*pow(xi, 3) + 144.0*pow(J1_, 3)*pow(J2_, 2)*k_s*pow(xi, 2) - 24.0*pow(J1_, 3)*pow(J2_, 2)*k_s + 240.0*pow(J1_, 3)*pow(J3_, 2)*k_s*pow(xi, 3) - 324.0*pow(J1_, 3)*pow(J3_, 2)*k_s*pow(xi, 2) - 72.0*pow(J1_, 3)*pow(J3_, 2)*k_s*xi + 54.0*pow(J1_, 3)*pow(J3_, 2)*k_s - 1440.0*pow(J1_, 3)*pow(k_s, 2)*pow(xi, 3) + 2160.0*pow(J1_, 3)*pow(k_s, 2)*pow(xi, 2) - 360.0*pow(J1_, 3)*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

d3Nt1_dxi: 
 (60*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*pow(xi, 2) - 24*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*xi - 6*pow(J1, 3)*pow(J2, 2)*pow(J3, 2) - 360.0*pow(J1, 3)*pow(J2, 2)*k_s*pow(xi, 2) + 288.0*pow(J1, 3)*pow(J2, 2)*k_s*xi + 720.0*pow(J1, 3)*pow(J3, 2)*k_s*pow(xi, 2) - 648.0*pow(J1, 3)*pow(J3, 2)*k_s*xi - 72.0*pow(J1, 3)*pow(J3, 2)*k_s - 4320.0*pow(J1, 3)*pow(k_s, 2)*pow(xi, 2) + 4320.0*pow(J1, 3)*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d4Nt1_dxi:
 (120*pow(J1, 3)*pow(J2, 2)*pow(J3, 2)*xi - 24*pow(J1, 3)*pow(J2, 2)*pow(J3, 2) - 720.0*pow(J1, 3)*pow(J2, 2)*k_s*xi + 288.0*pow(J1, 3)*pow(J2, 2)*k_s + 1440.0*pow(J1, 3)*pow(J3, 2)*k_s*xi - 648.0*pow(J1, 3)*pow(J3, 2)*k_s - 8640.0*pow(J1, 3)*pow(k_s, 2)*xi + 4320.0*pow(J1, 3)*pow(k_s, 2))/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))
'''

eq_n10 = sp.Eq(    u1, 0.0)
eq_n2  = sp.Eq(    v1, 0.0)
eq_n3  = sp.Eq(theta1, 0.0)
eq_n4  = sp.Eq(    u2, 1.0)
eq_n5  = sp.Eq(    v2, 0.0)
eq_n6  = sp.Eq(theta2, 0.0)
eq_n7  = sp.Eq(    u3, 0.0)
eq_n8  = sp.Eq(    v3, 0.0)
eq_n9  = sp.Eq(theta3, 0.0)
eq_n1  = sp.Eq(Nu2, poly)
solution_Nu2 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nu2, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nu2: \n", solution_Nu2[Nu2])
# print("dNu2_dxi: \n", sp.ccode(sp.diff(solution_Nu2[Nu2], xi)))
# print("d2Nu2_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nu2[Nu2], xi), xi)))
# print("d3Nu2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nu2[Nu2], xi), xi), xi)))
print("\n d4Nu2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nu2[Nu2], xi), xi), xi), xi)))

'''
Nu2:

 (-J1_**2*J2_**3*J3_**2*k02_*xi**5 - J1_**2*J2_**3*J3_**2*k02_*xi**4 + J1_**2*J2_**3*J3_**2*k02_*xi**3 + J1_**2*J2_**3*J3_**2*k02_*xi**2 + 6.0*J1_**2*J2_**3*k02_*k_s*xi**5 + 12.0*J1_**2*J2_**3*k02_*k_s*xi**4 - 12.0*J1_**2*J2_**3*k02_*k_s*xi**2 - 6.0*J1_**2*J2_**3*k02_*k_s*xi - 12.0*J2_**3*J3_**2*k02_*k_s*xi**5 - 27.0*J2_**3*J3_**2*k02_*k_s*xi**4 + 12.0*J2_**3*J3_**2*k02_*k_s*xi**3 + 27.0*J2_**3*J3_**2*k02_*k_s*xi**2 + 72.0*J2_**3*k02_*k_s**2*xi**5 + 180.0*J2_**3*k02_*k_s**2*xi**4 - 180.0*J2_**3*k02_*k_s**2*xi**2 - 72.0*J2_**3*k02_*k_s**2*xi)/(4.0*J1_**2*J2_**2*J3_**2 - 48.0*J1_**2*J2_**2*k_s + 78.0*J1_**2*J3_**2*k_s - 648.0*J1_**2*k_s**2 + 78.0*J2_**2*J3_**2*k_s - 648.0*J2_**2*k_s**2 + 1296.0*J3_**2*k_s**2 - 8640.0*k_s**3)

dNu2_dxi:
 (-5*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*pow(xi, 4) - 4*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*pow(xi, 3) + 3*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*pow(xi, 2) + 2*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*xi + 30.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s*pow(xi, 4) + 48.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s*pow(xi, 3) - 24.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s*xi - 6.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s - 60.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*pow(xi, 4) - 108.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*pow(xi, 3) + 36.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*pow(xi, 2) + 54.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*xi + 360.0*pow(J2_, 3)*k02_*pow(k_s, 2)*pow(xi, 4) + 720.0*pow(J2_, 3)*k02_*pow(k_s, 2)*pow(xi, 3) - 360.0*pow(J2_, 3)*k02_*pow(k_s, 2)*xi - 72.0*pow(J2_, 3)*k02_*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d2Nu2_dxi: 
 (-20*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*pow(xi, 3) - 12*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*pow(xi, 2) + 6*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_*xi + 2*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*k02_ + 120.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s*pow(xi, 3) + 144.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s*pow(xi, 2) - 24.0*pow(J1_, 2)*pow(J2_, 3)*k02_*k_s - 240.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*pow(xi, 3) - 324.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*pow(xi, 2) + 72.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s*xi + 54.0*pow(J2_, 3)*pow(J3_, 2)*k02_*k_s + 1440.0*pow(J2_, 3)*k02_*pow(k_s, 2)*pow(xi, 3) + 2160.0*pow(J2_, 3)*k02_*pow(k_s, 2)*pow(xi, 2) - 360.0*pow(J2_, 3)*k02_*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

d3Nu2_dxi: 
 (-60*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*k02*pow(xi, 2) - 24*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*k02*xi + 6*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*k02 + 360.0*pow(J1, 2)*pow(J2, 3)*k02*k_s*pow(xi, 2) + 288.0*pow(J1, 2)*pow(J2, 3)*k02*k_s*xi - 720.0*pow(J2, 3)*pow(J3, 2)*k02*k_s*pow(xi, 2) - 648.0*pow(J2, 3)*pow(J3, 2)*k02*k_s*xi + 72.0*pow(J2, 3)*pow(J3, 2)*k02*k_s + 4320.0*pow(J2, 3)*k02*pow(k_s, 2)*pow(xi, 2) + 4320.0*pow(J2, 3)*k02*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d4Nu2_dxi:
 
  (-120*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*k02*xi - 24*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*k02 + 720.0*pow(J1, 2)*pow(J2, 3)*k02*k_s*xi + 288.0*pow(J1, 2)*pow(J2, 3)*k02*k_s - 1440.0*pow(J2, 3)*pow(J3, 2)*k02*k_s*xi - 648.0*pow(J2, 3)*pow(J3, 2)*k02*k_s + 8640.0*pow(J2, 3)*k02*pow(k_s, 2)*xi + 4320.0*pow(J2, 3)*k02*pow(k_s, 2))/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

'''

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 1.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nv2, poly)
# solution_Nv2 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nv2, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nv2: \n", solution_Nv2[Nv2])
# print("dNv2_dxi: \n", sp.ccode(sp.diff(solution_Nv2[Nv2], xi)))
# print("d2Nv2_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nv2[Nv2], xi), xi)))
# print("d3Nv2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nv2[Nv2], xi), xi), xi)))
# print("\n d4Nv2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nv2[Nv2], xi), xi), xi), xi)))

'''
Nv2:

(-6.0*J1_**2*J2_**2*J3_**2*xi**5 - 4.0*J1_**2*J2_**2*J3_**2*xi**4 + 10.0*J1_**2*J2_**2*J3_**2*xi**3 + 8.0*J1_**2*J2_**2*J3_**2*xi**2 + 12.0*J1_**2*J2_**2*k_s*xi**5 + 48.0*J1_**2*J2_**2*k_s*xi**4 - 96.0*J1_**2*J2_**2*k_s*xi**2 - 60.0*J1_**2*J2_**2*k_s*xi - 18.0*J1_**2*J3_**2*k_s*xi**5 + 21.0*J1_**2*J3_**2*k_s*xi**4 + 96.0*J1_**2*J3_**2*k_s*xi**3 + 57.0*J1_**2*J3_**2*k_s*xi**2 - 72.0*J1_**2*k_s**2*xi**5 + 180.0*J1_**2*k_s**2*xi**4 - 828.0*J1_**2*k_s**2*xi**2 - 576.0*J1_**2*k_s**2*xi - 66.0*J2_**2*J3_**2*k_s*xi**5 - 129.0*J2_**2*J3_**2*k_s*xi**4 + 144.0*J2_**2*J3_**2*k_s*xi**3 + 207.0*J2_**2*J3_**2*k_s*xi**2 + 216.0*J2_**2*k_s**2*xi**5 + 540.0*J2_**2*k_s**2*xi**4 - 1188.0*J2_**2*k_s**2*xi**2 - 864.0*J2_**2*k_s**2*xi - 144.0*J3_**2*k_s**2*xi**5 + 1440.0*J3_**2*k_s**2*xi**3 + 1296.0*J3_**2*k_s**2*xi**2 - 8640.0*k_s**3*xi**2 - 8640.0*k_s**3*xi)/(8.0*J1_**2*J2_**2*J3_**2 - 96.0*J1_**2*J2_**2*k_s + 156.0*J1_**2*J3_**2*k_s - 1296.0*J1_**2*k_s**2 + 156.0*J2_**2*J3_**2*k_s - 1296.0*J2_**2*k_s**2 + 2592.0*J3_**2*k_s**2 - 17280.0*k_s**3)

dNv2_dxi: 
 (-30.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 4) - 16.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) + 30.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) + 16.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*xi + 60.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 4) + 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 3) - 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s*xi - 60.0*pow(J1_, 2)*pow(J2_, 2)*k_s - 90.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) + 84.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 288.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 114.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi - 360.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 
4) + 720.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) - 1656.0*pow(J1_, 2)*pow(k_s, 2)*xi - 576.0*pow(J1_, 2)*pow(k_s, 2) - 330.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) - 516.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 432.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 414.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi + 1080.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 4) + 2160.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) - 2376.0*pow(J2_, 2)*pow(k_s, 2)*xi - 864.0*pow(J2_, 2)*pow(k_s, 2) - 720.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 4) + 4320.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 2) + 2592.0*pow(J3_, 2)*pow(k_s, 2)*xi - 17280.0*pow(k_s, 3)*xi - 8640.0*pow(k_s, 3))/(8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 156.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J1_, 2)*pow(k_s, 2) + 156.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J2_, 2)*pow(k_s, 2) + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

d2Nv2_dxi: 
 (-120.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 48.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) + 60.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*xi + 16.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) + 240.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 3) + 576.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 2) - 192.0*pow(J1_, 2)*pow(J2_, 2)*k_s - 360.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 252.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 576.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi + 114.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 1440.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) + 2160.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 2) - 1656.0*pow(J1_, 2)*pow(k_s, 2) - 1320.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) - 1548.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 
2) + 864.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi + 414.0*pow(J2_, 2)*pow(J3_, 2)*k_s + 4320.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) + 6480.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 2) - 2376.0*pow(J2_, 2)*pow(k_s, 2) - 2880.0*pow(J3_, 2)*pow(k_s, 2)*pow(xi, 3) + 8640.0*pow(J3_, 2)*pow(k_s, 2)*xi + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))/(8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 156.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J1_, 2)*pow(k_s, 2) + 156.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 1296.0*pow(J2_, 
2)*pow(k_s, 2) + 2592.0*pow(J3_, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

d3Nv2_dxi: 
 (-360.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*pow(xi, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*xi + 60.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) + 720.0*pow(J1, 2)*pow(J2, 2)*k_s*pow(xi, 2) + 1152.0*pow(J1, 2)*pow(J2, 2)*k_s*xi - 1080.0*pow(J1, 2)*pow(J3, 2)*k_s*pow(xi, 2) + 504.0*pow(J1, 2)*pow(J3, 2)*k_s*xi + 576.0*pow(J1, 2)*pow(J3, 2)*k_s - 4320.0*pow(J1, 2)*pow(k_s, 2)*pow(xi, 2) + 4320.0*pow(J1, 2)*pow(k_s, 2)*xi - 3960.0*pow(J2, 2)*pow(J3, 2)*k_s*pow(xi, 2) - 3096.0*pow(J2, 2)*pow(J3, 2)*k_s*xi + 864.0*pow(J2, 2)*pow(J3, 2)*k_s + 12960.0*pow(J2, 2)*pow(k_s, 2)*pow(xi, 2) + 12960.0*pow(J2, 2)*pow(k_s, 2)*xi - 8640.0*pow(J3, 2)*pow(k_s, 2)*pow(xi, 2) + 8640.0*pow(J3, 2)*pow(k_s, 2))/(8.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*k_s + 156.0*pow(J1, 2)*pow(J3, 2)*k_s - 1296.0*pow(J1, 2)*pow(k_s, 2) + 156.0*pow(J2, 2)*pow(J3, 
2)*k_s - 1296.0*pow(J2, 2)*pow(k_s, 2) + 2592.0*pow(J3, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))

 d4Nv2_dxi: 
 (-720.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*xi - 96.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) + 1440.0*pow(J1, 2)*pow(J2, 2)*k_s*xi + 1152.0*pow(J1, 2)*pow(J2, 2)*k_s - 2160.0*pow(J1, 2)*pow(J3, 2)*k_s*xi + 504.0*pow(J1, 2)*pow(J3, 2)*k_s - 8640.0*pow(J1, 2)*pow(k_s, 2)*xi + 4320.0*pow(J1, 2)*pow(k_s, 2) - 7920.0*pow(J2, 2)*pow(J3, 2)*k_s*xi - 3096.0*pow(J2, 2)*pow(J3, 2)*k_s + 25920.0*pow(J2, 2)*pow(k_s, 2)*xi + 12960.0*pow(J2, 2)*pow(k_s, 2) - 17280.0*pow(J3, 2)*pow(k_s, 2)*xi)/(8.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 96.0*pow(J1, 2)*pow(J2, 2)*k_s + 156.0*pow(J1, 2)*pow(J3, 2)*k_s - 1296.0*pow(J1, 2)*pow(k_s, 2) + 156.0*pow(J2, 2)*pow(J3, 2)*k_s - 1296.0*pow(J2, 2)*pow(k_s, 2) + 2592.0*pow(J3, 2)*pow(k_s, 2) - 17280.0*pow(k_s, 3))
'''

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 1.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nt2, poly)
# solution_Nt2 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nt2, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nt2: \n", solution_Nt2[Nt2])
# print("dNt2_dxi: \n", sp.ccode(sp.diff(solution_Nt2[Nt2], xi)))
# print("d2Nt2_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nt2[Nt2], xi), xi)))
# print("d3Nt2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nt2[Nt2], xi), xi), xi)))
# print("\n d4Nt2_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nt2[Nt2], xi), xi), xi), xi)))


'''
Nt2:

 (J1_**2*J2_**3*J3_**2*xi**5 + J1_**2*J2_**3*J3_**2*xi**4 - J1_**2*J2_**3*J3_**2*xi**3 - J1_**2*J2_**3*J3_**2*xi**2 - 6.0*J1_**2*J2_**3*k_s*xi**5 - 12.0*J1_**2*J2_**3*k_s*xi**4 + 12.0*J1_**2*J2_**3*k_s*xi**2 + 6.0*J1_**2*J2_**3*k_s*xi + 12.0*J2_**3*J3_**2*k_s*xi**5 + 27.0*J2_**3*J3_**2*k_s*xi**4 - 12.0*J2_**3*J3_**2*k_s*xi**3 - 27.0*J2_**3*J3_**2*k_s*xi**2 - 72.0*J2_**3*k_s**2*xi**5 - 180.0*J2_**3*k_s**2*xi**4 + 180.0*J2_**3*k_s**2*xi**2 + 72.0*J2_**3*k_s**2*xi)/(4.0*J1_**2*J2_**2*J3_**2 - 48.0*J1_**2*J2_**2*k_s + 78.0*J1_**2*J3_**2*k_s - 648.0*J1_**2*k_s**2 + 78.0*J2_**2*J3_**2*k_s - 648.0*J2_**2*k_s**2 + 1296.0*J3_**2*k_s**2 - 8640.0*k_s**3)

 dNt2_dxi: 
 (5*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*pow(xi, 4) + 4*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*pow(xi, 3) - 3*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*pow(xi, 2) - 2*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*xi - 
30.0*pow(J1_, 2)*pow(J2_, 3)*k_s*pow(xi, 4) - 48.0*pow(J1_, 2)*pow(J2_, 3)*k_s*pow(xi, 3) + 24.0*pow(J1_, 2)*pow(J2_, 3)*k_s*xi + 6.0*pow(J1_, 2)*pow(J2_, 3)*k_s + 60.0*pow(J2_, 3)*pow(J3_, 2)*k_s*pow(xi, 4) + 108.0*pow(J2_, 3)*pow(J3_, 2)*k_s*pow(xi, 3) - 36.0*pow(J2_, 3)*pow(J3_, 2)*k_s*pow(xi, 2) - 54.0*pow(J2_, 3)*pow(J3_, 2)*k_s*xi - 360.0*pow(J2_, 3)*pow(k_s, 2)*pow(xi, 4) - 720.0*pow(J2_, 3)*pow(k_s, 2)*pow(xi, 3) + 360.0*pow(J2_, 3)*pow(k_s, 2)*xi + 72.0*pow(J2_, 3)*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

d2Nt2_dxi: 
 (20*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*pow(xi, 3) + 12*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*pow(xi, 2) - 6*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2)*xi - 2*pow(J1_, 2)*pow(J2_, 3)*pow(J3_, 2) - 120.0*pow(J1_, 2)*pow(J2_, 3)*k_s*pow(xi, 3) - 144.0*pow(J1_, 2)*pow(J2_, 3)*k_s*pow(xi, 2) + 24.0*pow(J1_, 2)*pow(J2_, 3)*k_s + 240.0*pow(J2_, 3)*pow(J3_, 2)*k_s*pow(xi, 3) + 324.0*pow(J2_, 3)*pow(J3_, 2)*k_s*pow(xi, 2) - 72.0*pow(J2_, 3)*pow(J3_, 2)*k_s*xi - 54.0*pow(J2_, 3)*pow(J3_, 2)*k_s - 1440.0*pow(J2_, 3)*pow(k_s, 2)*pow(xi, 3) - 2160.0*pow(J2_, 3)*pow(k_s, 2)*pow(xi, 2) + 360.0*pow(J2_, 3)*pow(k_s, 2))/(4.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 48.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 78.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J1_, 2)*pow(k_s, 2) + 78.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 648.0*pow(J2_, 2)*pow(k_s, 2) + 1296.0*pow(J3_, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d3Nt2_dxi: 
 (60*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*pow(xi, 2) + 24*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*xi - 6*pow(J1, 2)*pow(J2, 3)*pow(J3, 2) - 360.0*pow(J1, 2)*pow(J2, 3)*k_s*pow(xi, 2) - 288.0*pow(J1, 2)*pow(J2, 3)*k_s*xi + 720.0*pow(J2, 3)*pow(J3, 2)*k_s*pow(xi, 2) + 648.0*pow(J2, 3)*pow(J3, 2)*k_s*xi - 72.0*pow(J2, 3)*pow(J3, 2)*k_s - 4320.0*pow(J2, 3)*pow(k_s, 2)*pow(xi, 2) - 4320.0*pow(J2, 3)*pow(k_s, 2)*xi)/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) + 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))

 d4Nt2_dxi: 
 (120*pow(J1, 2)*pow(J2, 3)*pow(J3, 2)*xi + 24*pow(J1, 2)*pow(J2, 3)*pow(J3, 2) - 720.0*pow(J1, 2)*pow(J2, 3)*k_s*xi - 288.0*pow(J1, 2)*pow(J2, 3)*k_s + 1440.0*pow(J2, 3)*pow(J3, 2)*k_s*xi + 648.0*pow(J2, 3)*pow(J3, 2)*k_s - 8640.0*pow(J2, 3)*pow(k_s, 2)*xi - 4320.0*pow(J2, 3)*pow(k_s, 2))/(4.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 48.0*pow(J1, 2)*pow(J2, 2)*k_s + 78.0*pow(J1, 2)*pow(J3, 2)*k_s - 648.0*pow(J1, 2)*pow(k_s, 2) 
+ 78.0*pow(J2, 2)*pow(J3, 2)*k_s - 648.0*pow(J2, 2)*pow(k_s, 2) + 1296.0*pow(J3, 2)*pow(k_s, 2) - 8640.0*pow(k_s, 3))
'''

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 1.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nu3, poly)
# solution_Nu3 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nu3, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nu3: \n", sp.ccode(solution_Nu3[Nu3]))
# print("dNu3_dxi: \n", sp.ccode(sp.diff(solution_Nu3[Nu3], xi)))
# print("d2Nu3_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nu3[Nu3], xi), xi)))
# print("d3Nu3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nu3[Nu3], xi), xi), xi)))
# print("\n d4Nu3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nu3[Nu3], xi), xi), xi), xi)))


'''
Nu3: 

 (-2.0*J1_**2*J2_**2*J3_**3*k03_*xi**5 + 4.0*J1_**2*J2_**2*J3_**3*k03_*xi**3 - 2.0*J1_**2*J2_**2*J3_**3*k03_*xi - 15.0*J1_**2*J3_**3*k03_*k_s*xi**5 + 24.0*J1_**2*J3_**3*k03_*k_s*xi**4 + 54.0*J1_**2*J3_**3*k03_*k_s*xi**3 - 24.0*J1_**2*J3_**3*k03_*k_s*xi**2 - 39.0*J1_**2*J3_**3*k03_*k_s*xi - 15.0*J2_**2*J3_**3*k03_*k_s*xi**5 - 24.0*J2_**2*J3_**3*k03_*k_s*xi**4 + 54.0*J2_**2*J3_**3*k03_*k_s*xi**3 + 24.0*J2_**2*J3_**3*k03_*k_s*xi**2 - 39.0*J2_**2*J3_**3*k03_*k_s*xi - 72.0*J3_**3*k03_*k_s**2*xi**5 + 720.0*J3_**3*k03_*k_s**2*xi**3 - 648.0*J3_**3*k03_*k_s**2*xi)/(2.0*J1_**2*J2_**2*J3_**2 - 24.0*J1_**2*J2_**2*k_s + 
39.0*J1_**2*J3_**2*k_s - 324.0*J1_**2*k_s**2 + 39.0*J2_**2*J3_**2*k_s - 324.0*J2_**2*k_s**2 + 648.0*J3_**2*k_s**2 - 4320.0*k_s**3)

dNu3_dxi: 
 (-10.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03*pow(xi, 4) + 12.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03*pow(xi, 2) - 2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03 - 75.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*pow(xi, 4) + 96.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*pow(xi, 3) + 162.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*pow(xi, 2) - 48.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*xi - 39.0*pow(J1, 2)*pow(J3, 3)*k03*k_s - 75.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*pow(xi, 4) - 96.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*pow(xi, 3) + 162.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*pow(xi, 2) + 48.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*xi - 39.0*pow(J2, 2)*pow(J3, 3)*k03*k_s - 360.0*pow(J3, 3)*k03*pow(k_s, 2)*pow(xi, 4) + 2160.0*pow(J3, 3)*k03*pow(k_s, 2)*pow(xi, 2) - 648.0*pow(J3, 3)*k03*pow(k_s, 2))/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

d2Nu3_dxi: 
 (-40.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 3)*k03_*pow(xi, 3) + 24.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 3)*k03_*xi - 300.0*pow(J1_, 2)*pow(J3_, 3)*k03_*k_s*pow(xi, 3) + 288.0*pow(J1_, 2)*pow(J3_, 3)*k03_*k_s*pow(xi, 2) + 324.0*pow(J1_, 2)*pow(J3_, 3)*k03_*k_s*xi - 48.0*pow(J1_, 2)*pow(J3_, 3)*k03_*k_s - 300.0*pow(J2_, 2)*pow(J3_, 3)*k03_*k_s*pow(xi, 3) - 288.0*pow(J2_, 2)*pow(J3_, 3)*k03_*k_s*pow(xi, 2) + 324.0*pow(J2_, 2)*pow(J3_, 3)*k03_*k_s*xi + 48.0*pow(J2_, 2)*pow(J3_, 3)*k03_*k_s - 1440.0*pow(J3_, 3)*k03_*pow(k_s, 2)*pow(xi, 3) + 4320.0*pow(J3_, 3)*k03_*pow(k_s, 2)*xi)/(2.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 24.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 39.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J1_, 2)*pow(k_s, 2) + 39.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J2_, 2)*pow(k_s, 2) + 648.0*pow(J3_, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

d3Nu3_dxi: 
 (-120.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03*pow(xi, 2) + 24.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03 - 900.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*pow(xi, 2) + 576.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*xi + 324.0*pow(J1, 2)*pow(J3, 3)*k03*k_s - 900.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*pow(xi, 2) - 576.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*xi + 324.0*pow(J2, 2)*pow(J3, 3)*k03*k_s - 4320.0*pow(J3, 3)*k03*pow(k_s, 2)*pow(xi, 2) + 4320.0*pow(J3, 3)*k03*pow(k_s, 2))/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

 d4Nu3_dxi:
 (-240.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*k03*xi - 1800.0*pow(J1, 2)*pow(J3, 3)*k03*k_s*xi + 576.0*pow(J1, 2)*pow(J3, 3)*k03*k_s - 1800.0*pow(J2, 2)*pow(J3, 3)*k03*k_s*xi - 576.0*pow(J2, 2)*pow(J3, 3)*k03*k_s - 8640.0*pow(J3, 3)*k03*pow(k_s, 2)*xi)/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 
2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))
'''


# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 1.0)
# eq_n9  = sp.Eq(theta3, 0.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nv3, poly)
# solution_Nv3 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nv3, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nv3: \n", solution_Nv3[Nv3])
# print("dNv3_dxi: \n", sp.ccode(sp.diff(solution_Nv3[Nv3], xi)))
# print("d2Nv3_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nv3[Nv3], xi), xi)))
# print("d3Nv3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nv3[Nv3], xi), xi), xi)))
# print("\n d4Nv3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nv3[Nv3], xi), xi), xi), xi)))

'''
Nv3: 
 (2.0*J1_**2*J2_**2*J3_**2*xi**4 - 4.0*J1_**2*J2_**2*J3_**2*xi**2 + 2.0*J1_**2*J2_**2*J3_**2 - 24.0*J1_**2*J2_**2*k_s*xi**4 + 48.0*J1_**2*J2_**2*k_s*xi**2 - 24.0*J1_**2*J2_**2*k_s - 12.0*J1_**2*J3_**2*k_s*xi**5 + 27.0*J1_**2*J3_**2*k_s*xi**4 + 12.0*J1_**2*J3_**2*k_s*xi**3 - 66.0*J1_**2*J3_**2*k_s*xi**2 + 39.0*J1_**2*J3_**2*k_s + 72.0*J1_**2*k_s**2*xi**5 - 180.0*J1_**2*k_s**2*xi**4 + 504.0*J1_**2*k_s**2*xi**2 - 72.0*J1_**2*k_s**2*xi - 324.0*J1_**2*k_s**2 + 12.0*J2_**2*J3_**2*k_s*xi**5 + 27.0*J2_**2*J3_**2*k_s*xi**4 - 12.0*J2_**2*J3_**2*k_s*xi**3 - 66.0*J2_**2*J3_**2*k_s*xi**2 + 39.0*J2_**2*J3_**2*k_s - 72.0*J2_**2*k_s**2*xi**5 - 180.0*J2_**2*k_s**2*xi**4 + 504.0*J2_**2*k_s**2*xi**2 + 72.0*J2_**2*k_s**2*xi - 324.0*J2_**2*k_s**2 - 648.0*J3_**2*k_s**2*xi**2 + 648.0*J3_**2*k_s**2 + 4320.0*k_s**3*xi**2 - 4320.0*k_s**3)/(2.0*J1_**2*J2_**2*J3_**2 - 24.0*J1_**2*J2_**2*k_s + 39.0*J1_**2*J3_**2*k_s - 324.0*J1_**2*k_s**2 + 39.0*J2_**2*J3_**2*k_s - 324.0*J2_**2*k_s**2 + 648.0*J3_**2*k_s**2 - 4320.0*k_s**3)

 dNv3_dxi: 
 (8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 3) - 8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*xi - 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 3) + 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s*xi - 60.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) + 108.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 36.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) - 132.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi + 360.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 4) - 720.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) + 1008.0*pow(J1_, 2)*pow(k_s, 2)*xi - 72.0*pow(J1_, 2)*pow(k_s, 2) + 60.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 4) + 108.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) - 36.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) - 132.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi - 360.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 4) - 720.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) + 1008.0*pow(J2_, 2)*pow(k_s, 2)*xi + 72.0*pow(J2_, 2)*pow(k_s, 2) - 1296.0*pow(J3_, 2)*pow(k_s, 2)*xi + 8640.0*pow(k_s, 3)*xi)/(2.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 24.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 39.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J1_, 2)*pow(k_s, 2) + 39.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J2_, 2)*pow(k_s, 2) + 648.0*pow(J3_, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

 d2Nv3_dxi: 
 (24.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2)*pow(xi, 2) - 8.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 288.0*pow(J1_, 2)*pow(J2_, 2)*k_s*pow(xi, 2) + 96.0*pow(J1_, 2)*pow(J2_, 2)*k_s 
- 240.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 324.0*pow(J1_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) + 72.0*pow(J1_, 2)*pow(J3_, 2)*k_s*xi - 132.0*pow(J1_, 2)*pow(J3_, 2)*k_s + 1440.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 3) - 2160.0*pow(J1_, 2)*pow(k_s, 2)*pow(xi, 2) + 1008.0*pow(J1_, 2)*pow(k_s, 2) + 240.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 3) + 324.0*pow(J2_, 2)*pow(J3_, 2)*k_s*pow(xi, 2) - 72.0*pow(J2_, 2)*pow(J3_, 2)*k_s*xi - 132.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 1440.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 3) - 2160.0*pow(J2_, 2)*pow(k_s, 2)*pow(xi, 2) + 1008.0*pow(J2_, 2)*pow(k_s, 2) - 1296.0*pow(J3_, 2)*pow(k_s, 2) + 8640.0*pow(k_s, 3))/(2.0*pow(J1_, 2)*pow(J2_, 2)*pow(J3_, 2) - 24.0*pow(J1_, 2)*pow(J2_, 2)*k_s + 39.0*pow(J1_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J1_, 2)*pow(k_s, 2) + 39.0*pow(J2_, 2)*pow(J3_, 2)*k_s - 324.0*pow(J2_, 2)*pow(k_s, 2) + 648.0*pow(J3_, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

d3Nv3_dxi: 
 (48.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2)*xi - 576.0*pow(J1, 2)*pow(J2, 2)*k_s*xi - 720.0*pow(J1, 2)*pow(J3, 2)*k_s*pow(xi, 2) + 648.0*pow(J1, 2)*pow(J3, 2)*k_s*xi + 72.0*pow(J1, 2)*pow(J3, 2)*k_s + 4320.0*pow(J1, 2)*pow(k_s, 2)*pow(xi, 2) - 4320.0*pow(J1, 2)*pow(k_s, 2)*xi + 720.0*pow(J2, 2)*pow(J3, 2)*k_s*pow(xi, 2) + 648.0*pow(J2, 2)*pow(J3, 2)*k_s*xi - 72.0*pow(J2, 2)*pow(J3, 2)*k_s - 4320.0*pow(J2, 2)*pow(k_s, 2)*pow(xi, 2) - 4320.0*pow(J2, 2)*pow(k_s, 2)*xi)/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

 d4Nv3_dxi:
 (48.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 576.0*pow(J1, 2)*pow(J2, 2)*k_s - 1440.0*pow(J1, 2)*pow(J3, 2)*k_s*xi + 648.0*pow(J1, 2)*pow(J3, 2)*k_s + 8640.0*pow(J1, 2)*pow(k_s, 2)*xi - 4320.0*pow(J1, 2)*pow(k_s, 2) + 1440.0*pow(J2, 2)*pow(J3, 2)*k_s*xi + 648.0*pow(J2, 2)*pow(J3, 2)*k_s - 8640.0*pow(J2, 2)*pow(k_s, 2)*xi - 4320.0*pow(J2, 2)*pow(k_s, 2))/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 
39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 
2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))
'''

# eq_n2  = sp.Eq(    v1, 0.0)
# eq_n3  = sp.Eq(theta1, 0.0)
# eq_n4  = sp.Eq(    u2, 0.0)
# eq_n5  = sp.Eq(    v2, 0.0)
# eq_n6  = sp.Eq(theta2, 0.0)
# eq_n7  = sp.Eq(    u3, 0.0)
# eq_n8  = sp.Eq(    v3, 0.0)
# eq_n9  = sp.Eq(theta3, 1.0)
# eq_n10 = sp.Eq(    u1, 0.0)
# eq_n1  = sp.Eq(Nt3, poly)
# solution_Nt3 = sp.solve((eq_n10, eq_n2, eq_n3, eq_n4, eq_n5, eq_n6, eq_n7, eq_n8, eq_n9, eq_n1), (Nt3, u1, v1, theta1, u2, v2, theta2, u3, v3, theta3))
# print("Nt3: \n", solution_Nt3[Nt3])
# print("dNt3_dxi: \n", sp.ccode(sp.diff(solution_Nt3[Nt3], xi)))
# print("d2Nt3_dxi: \n", sp.ccode(sp.diff(sp.diff(solution_Nt3[Nt3], xi), xi)))
# print("d3Nt3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(solution_Nt3[Nt3], xi), xi), xi)))
# print("\n d4Nt3_dxi: \n", sp.ccode(sp.diff(sp.diff(sp.diff(sp.diff(solution_Nt3[Nt3], xi), xi), xi), xi)))

'''
Nt3: 
 (2.0*J1_**2*J2_**2*J3_**3*xi**5 - 4.0*J1_**2*J2_**2*J3_**3*xi**3 + 2.0*J1_**2*J2_**2*J3_**3*xi + 15.0*J1_**2*J3_**3*k_s*xi**5 - 24.0*J1_**2*J3_**3*k_s*xi**4 - 54.0*J1_**2*J3_**3*k_s*xi**3 + 24.0*J1_**2*J3_**3*k_s*xi**2 + 39.0*J1_**2*J3_**3*k_s*xi + 15.0*J2_**2*J3_**3*k_s*xi**5 + 24.0*J2_**2*J3_**3*k_s*xi**4 - 54.0*J2_**2*J3_**3*k_s*xi**3 - 24.0*J2_**2*J3_**3*k_s*xi**2 + 39.0*J2_**2*J3_**3*k_s*xi + 72.0*J3_**3*k_s**2*xi**5 - 720.0*J3_**3*k_s**2*xi**3 + 648.0*J3_**3*k_s**2*xi)/(2.0*J1_**2*J2_**2*J3_**2 - 24.0*J1_**2*J2_**2*k_s + 39.0*J1_**2*J3_**2*k_s - 324.0*J1_**2*k_s**2 + 39.0*J2_**2*J3_**2*k_s - 324.0*J2_**2*k_s**2 + 648.0*J3_**2*k_s**2 - 4320.0*k_s**3)

 dNt3_dxi: 
(10.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*pow(xi, 4) - 12.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*pow(xi, 2) + 2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3) + 75.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 4) - 96.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 3) - 162.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 2) + 48.0*pow(J1, 2)*pow(J3, 3)*k_s*xi + 39.0*pow(J1, 2)*pow(J3, 3)*k_s + 75.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 4) + 96.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 3) - 162.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 2) - 48.0*pow(J2, 2)*pow(J3, 3)*k_s*xi + 39.0*pow(J2, 2)*pow(J3, 3)*k_s + 360.0*pow(J3, 3)*pow(k_s, 2)*pow(xi, 4) - 2160.0*pow(J3, 3)*pow(k_s, 2)*pow(xi, 2) + 648.0*pow(J3, 3)*pow(k_s, 2))/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

d2Nt3_dxi: 
 (40.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*pow(xi, 3) - 24.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*xi + 300.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 3) - 288.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 2) - 324.0*pow(J1, 2)*pow(J3, 3)*k_s*xi + 48.0*pow(J1, 2)*pow(J3, 3)*k_s + 300.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 3) + 288.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 2) - 324.0*pow(J2, 2)*pow(J3, 3)*k_s*xi - 48.0*pow(J2, 2)*pow(J3, 3)*k_s + 1440.0*pow(J3, 3)*pow(k_s, 2)*pow(xi, 3) - 4320.0*pow(J3, 3)*pow(k_s, 2)*xi)/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

 d3Nt3_dxi: 
 (120.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*pow(xi, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3) + 900.0*pow(J1, 2)*pow(J3, 3)*k_s*pow(xi, 2) - 576.0*pow(J1, 2)*pow(J3, 3)*k_s*xi - 324.0*pow(J1, 2)*pow(J3, 3)*k_s + 900.0*pow(J2, 2)*pow(J3, 3)*k_s*pow(xi, 2) + 576.0*pow(J2, 2)*pow(J3, 3)*k_s*xi - 324.0*pow(J2, 2)*pow(J3, 3)*k_s + 4320.0*pow(J3, 3)*pow(k_s, 2)*pow(xi, 2) - 4320.0*pow(J3, 3)*pow(k_s, 2))/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s 
- 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))

 d4Nt3_dxi:
 (240.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 3)*xi + 1800.0*pow(J1, 2)*pow(J3, 3)*k_s*xi - 576.0*pow(J1, 2)*pow(J3, 3)*k_s + 1800.0*pow(J2, 2)*pow(J3, 3)*k_s*xi + 576.0*pow(J2, 2)*pow(J3, 3)*k_s + 8640.0*pow(J3, 3)*pow(k_s, 2)*xi)/(2.0*pow(J1, 2)*pow(J2, 2)*pow(J3, 2) - 24.0*pow(J1, 2)*pow(J2, 2)*k_s + 39.0*pow(J1, 2)*pow(J3, 2)*k_s - 324.0*pow(J1, 2)*pow(k_s, 2) + 39.0*pow(J2, 2)*pow(J3, 2)*k_s - 324.0*pow(J2, 2)*pow(k_s, 2) + 648.0*pow(J3, 2)*pow(k_s, 2) - 4320.0*pow(k_s, 3))
'''
