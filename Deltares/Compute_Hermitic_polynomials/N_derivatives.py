import sympy as sp

# Define symbols
phi, xi, L = sp.symbols('phi xi L')

# Define N_i expressions
N = [
    (-40.0*phi**2 - 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
        (16.0*phi + 8.0 )/(32.0*phi + 8.0) * xi**2 + \
        (40.0*phi + 10.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
        (- 4.0 )/(32.0*phi + 8.0) * xi**4 + \
        (- 6.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5,

    (-L*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
        (L)/(32.0*phi + 8.0) * xi**2 + \
        (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
        (-L)/(32.0*phi + 8.0) * xi**4 + \
        (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5,

    1.0 + \
      (- 32.0*phi  - 16.0 )/(32.0*phi + 8.0) * xi**2 + \
      (8.0)/(32.0*phi + 8.0) * xi**4,

    (- 18.0*L*phi - 2.0*L ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
         (40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
         (- 4.0*L*phi  - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5,

    (40.0*phi**2 + 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
        (16.0*phi + 8.0)/(32.0*phi + 8.0) * xi**2 + \
        (- 40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 +\
        (- 4.0)/(32.0*phi + 8.0) * xi**4 +\
        (6.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5,

    (- L*phi ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
         (- L)/(32.0*phi + 8.0) * xi**2 + \
         (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
         (L)/(32.0*phi + 8.0) * xi**4 + \
         (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5
]

# # Compute derivatives of each N_i with respect to xi
# N_derivatives = [sp.diff(N_i, xi) for N_i in N]

# # Print the derivatives
# for i, N_derivative in enumerate(N_derivatives):
#     print(f"d(N[{i}])/d(xi) =", N_derivative)


'''
dN_dxi:

    d(N[0])/d(xi) = -30.0*xi**4/(160.0*phi**2 - 40.0*phi - 8.0) - 16.0*xi**3/(32.0*phi + 8.0) + 3*xi**2*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*xi*(16.0*phi + 8.0)/(32.0*phi + 8.0) + (-40.0*phi**2 - 10.0*phi)/(80.0*phi**2 - 20.0*phi - 4.0)
    d(N[1])/d(xi) = -L*phi/(80.0*phi**2 - 20.0*phi - 4.0) - 4*L*xi**3/(32.0*phi + 8.0) + 3*L*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) + 2*L*xi/(32.0*phi + 8.0) + 5*xi**4*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
    d(N[2])/d(xi) = 32.0*xi**3/(32.0*phi + 8.0) + 2*xi*(-32.0*phi - 16.0)/(32.0*phi + 8.0)
    d(N[3])/d(xi) = 5*xi**4*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 3*xi**2*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + (-18.0*L*phi - 2.0*L)/(80.0*phi**2 - 20.0*phi - 4.0)
    d(N[4])/d(xi) = 30.0*xi**4/(160.0*phi**2 - 40.0*phi - 8.0) - 16.0*xi**3/(32.0*phi + 8.0) + 3*xi**2*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*xi*(16.0*phi + 8.0)/(32.0*phi + 8.0) + (40.0*phi**2 + 10.0*phi)/(80.0*phi**2 - 20.0*phi - 4.0)
    d(N[5])/d(xi) = -L*phi/(80.0*phi**2 - 20.0*phi - 4.0) + 4*L*xi**3/(32.0*phi + 8.0) + 3*L*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 2*L*xi/(32.0*phi + 8.0) + 5*xi**4*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
'''

# # Compute derivatives of each N_i with respect to xi
# N_second_derivatives = [sp.diff(sp.diff(N_i, xi), xi) for N_i in N]

# # Print the derivatives
# for i, N_second_derivatives in enumerate(N_second_derivatives):
#     print(f"d2(N[{i}])/d(xi2) =", N_second_derivatives)


'''
d2N_dxi2:

d2(N[0])/d(xi2) = -120.0*xi**3/(160.0*phi**2 - 40.0*phi - 8.0) - 48.0*xi**2/(32.0*phi + 8.0) + 6*xi*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*(16.0*phi + 8.0)/(32.0*phi + 8.0)
d2(N[1])/d(xi2) = -12*L*xi**2/(32.0*phi + 8.0) + 6*L*xi/(160.0*phi**2 - 40.0*phi - 8.0) + 2*L/(32.0*phi + 8.0) + 20*xi**3*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
d2(N[2])/d(xi2) = 96.0*xi**2/(32.0*phi + 8.0) + 2*(-32.0*phi - 16.0)/(32.0*phi + 8.0)
d2(N[3])/d(xi2) = 20*xi**3*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 6*xi*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)
d2(N[4])/d(xi2) = 120.0*xi**3/(160.0*phi**2 - 40.0*phi - 8.0) - 48.0*xi**2/(32.0*phi + 8.0) + 6*xi*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*(16.0*phi + 8.0)/(32.0*phi + 8.0)
d2(N[5])/d(xi2) = 12*L*xi**2/(32.0*phi + 8.0) + 6*L*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 2*L/(32.0*phi + 8.0) + 20*xi**3*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
'''

# # Compute derivatives of each N_i with respect to xi
# N_third_derivatives = [sp.diff(sp.diff(sp.diff(N_i, xi), xi), xi) for N_i in N]

# # Print the derivatives
# for i, N_third_derivatives in enumerate(N_third_derivatives):
#     print(f"d3(N[{i}])/d(xi3) =", N_third_derivatives)


'''
d3N_dxi3:

d3(N[0])/d(xi3) = -360.0*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0*xi/(32.0*phi + 8.0) + 6*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0)
d3(N[1])/d(xi3) = -24*L*xi/(32.0*phi + 8.0) + 6*L/(160.0*phi**2 - 40.0*phi - 8.0) + 60*xi**2*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)   
d3(N[2])/d(xi3) = 192.0*xi/(32.0*phi + 8.0)
d3(N[3])/d(xi3) = 60*xi**2*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 6*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)       
d3(N[4])/d(xi3) = 360.0*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0*xi/(32.0*phi + 8.0) + 6*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0)
d3(N[5])/d(xi3) = 24*L*xi/(32.0*phi + 8.0) + 6*L/(160.0*phi**2 - 40.0*phi - 8.0) + 60*xi**2*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
'''

# Compute derivatives of each N_i with respect to xi
N_fourth_derivatives = [sp.diff(sp.diff(sp.diff(sp.diff(N_i, xi), xi), xi), xi) for N_i in N]

# Print the derivatives
for i, N_fourth_derivatives in enumerate(N_fourth_derivatives):
    print(f"d4(N[{i}])/d(xi4) =", N_fourth_derivatives)

'''
d4(N[0])/d(xi4) = -720.0*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0/(32.0*phi + 8.0)
d4(N[1])/d(xi4) = -24*L/(32.0*phi + 8.0) + 120*xi*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
d4(N[2])/d(xi4) = 192.0/(32.0*phi + 8.0)
d4(N[3])/d(xi4) = 120*xi*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)
d4(N[4])/d(xi4) = 720.0*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0/(32.0*phi + 8.0)
d4(N[5])/d(xi4) = 24*L/(32.0*phi + 8.0) + 120*xi*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
'''