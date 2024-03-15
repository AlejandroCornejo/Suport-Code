import timoshenko_element_2D_2N
import math
import numpy as np

"""  

Felippa and Oñate, "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability", Archives of Comp. Methods in Eng. (2021) 28:2021-2080

In this file we derive the main equations of 3-noded Timoshenko FE, using a 5th order hermitic polynomials (derived by A. Cornejo)

Notation:

- El        : Axial strain (du0 / dx)
- Gamma_xy  : Shear strain
- Kappa     : curvature (d Theta / dx)
- E         : Young modulus
- A         : Cross section area
- I         : Inertia
- nu        : Poisson ratio
- G         : Shear modulus
- Phi       : Shear slenderness
- J         : Jacobian to the isoparametric space

- Nu        : Axial displacement shape functions
- N         : Transverse displacement shape functions

--> 9 Dofs per element: [u01, v1, theta1,   u02, v2, theta2,   u03, v3, theta3]

"""

class TimoshenkoElement2D3N(timoshenko_element_2D_2N.TimoshenkoElement2D2N):
    # ------------------------------------------------------------------------------------------------
    def __init__(self, Node1, Node2, Node3, E, I, nu, A):
        self.Nodes = [Node1, Node2, Node3] # The intermediate is the midpoint node
        self.E      = E
        self.A      = A
        self.As     = A * 5.0 / 6.0 # NOTE depends on the cross section
        self.I      = I
        self.nu     = nu
        self.Length = math.sqrt((Node1.x - Node3.x)**2 + (Node1.y - Node3.y)**2)
        self.alpha  = self.GetAlphaAngle()
        self.G      = E / (2.0 * (1.0 + nu))
        self.Phi    = 12.0 * E * I / (self.G * self.As * self.Length**2)
        self.IntegrationOrder = 2
        self.J      = self.Length * 0.5
    # ------------------------------------------------------------------------------------------------
    def GetShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Precompute powers of xi
        xi_sq = xi ** 2
        xi_cube = xi ** 3
        xi_quad = xi ** 4
        xi_quint = xi ** 5

        # Common denominators
        denom1 = 80.0 * phi ** 2 - 20.0 * phi - 4.0
        denom2 = 32.0 * phi + 8.0
        denom3 = 160.0 * phi ** 2 - 40.0 * phi - 8.0

        # Shape function expressions
        shape_functions_values = np.zeros(6)

        shape_functions_values[0] = (-40.0 * phi ** 2 - 10.0 * phi) / denom1 * xi + \
                                    (16.0 * phi + 8.0) / denom2 * xi_sq + \
                                    (40.0 * phi + 10.0) / denom3 * xi_cube + \
                                    (-4.0) / denom2 * xi_quad + \
                                    (-6.0) / denom3 * xi_quint

        shape_functions_values[1] = (-L * phi) / denom1 * xi + \
                                    L / denom2 * xi_sq + \
                                    L / denom3 * xi_cube + \
                                    (-L) / denom2 * xi_quad + \
                                    (2.0 * L * phi - L) / denom3 * xi_quint

        shape_functions_values[2] = 1.0 + \
                                    (-32.0 * phi - 16.0) / denom2 * xi_sq + \
                                    8.0 / denom2 * xi_quad

        shape_functions_values[3] = (-18.0 * L * phi - 2.0 * L) / denom1 * xi + \
                                    (40.0 * L * phi + 8.0 * L) / denom3 * xi_cube + \
                                    (-4.0 * L * phi - 4.0 * L) / denom3 * xi_quint

        shape_functions_values[4] = (40.0 * phi ** 2 + 10.0 * phi) / denom1 * xi + \
                                    (16.0 * phi + 8.0) / denom2 * xi_sq + \
                                    (-40.0 * phi - 10.0) / denom3 * xi_cube + \
                                    (-4.0) / denom2 * xi_quad + \
                                    6.0 / denom3 * xi_quint

        shape_functions_values[5] = (-L * phi) / denom1 * xi + \
                                    (-L) / denom2 * xi_sq + \
                                    L / denom3 * xi_cube + \
                                    L / denom2 * xi_quad + \
                                    (2.0 * L * phi - L) / denom3 * xi_quint
        return shape_functions_values

    # ------------------------------------------------------------------------------------------------
    def GetFirstDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Precompute powers of xi
        xi_sq = xi ** 2
        xi_cube = xi ** 3
        xi_quad = xi ** 4

        # Common denominators
        denom1 = 80.0 * phi ** 2 - 20.0 * phi - 4.0
        denom2 = 32.0 * phi + 8.0
        denom3 = 160.0 * phi ** 2 - 40.0 * phi - 8.0

        # Shape function derivatives expressions
        shape_functions_derivatives_values = np.zeros(6)

        shape_functions_derivatives_values[0] = -30.0 * xi_quad / denom3 - \
                                                16.0 * xi_cube / denom2 + \
                                                3 * xi_sq * (40.0 * phi + 10.0) / denom3 + \
                                                2 * xi * (16.0 * phi + 8.0) / denom2 - \
                                                (-40.0 * phi ** 2 - 10.0 * phi) / denom1

        shape_functions_derivatives_values[1] = -L * phi / denom1 - \
                                                4 * L * xi_cube / denom2 + \
                                                3 * L * xi_sq / denom3 + \
                                                2 * L * xi / denom2 + \
                                                5 * xi_quad * (2.0 * L * phi - L) / denom3

        shape_functions_derivatives_values[2] = 32.0 * xi_cube / denom2 + \
                                                2 * xi * (-32.0 * phi - 16.0) / denom2

        shape_functions_derivatives_values[3] = 5 * xi_quad * (-4.0 * L * phi - 4.0 * L) / denom3 + \
                                                3 * xi_sq * (40.0 * L * phi + 8.0 * L) / denom3 - \
                                                (-18.0 * L * phi - 2.0 * L) / denom1

        shape_functions_derivatives_values[4] = 30.0 * xi_quad / denom3 - \
                                                16.0 * xi_cube / denom2 + \
                                                3 * xi_sq * (-40.0 * phi - 10.0) / denom3 + \
                                                2 * xi * (16.0 * phi + 8.0) / denom2 + \
                                                (40.0 * phi ** 2 + 10.0 * phi) / denom1

        shape_functions_derivatives_values[5] = -L * phi / denom1 + \
                                                4 * L * xi_cube / denom2 + \
                                                3 * L * xi_sq / denom3 - \
                                                2 * L * xi / denom2 + \
                                                5 * xi_quad * (2.0 * L * phi - L) / denom3
        return shape_functions_derivatives_values
    # ------------------------------------------------------------------------------------------------
    def GetSecondDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Precompute powers of xi
        xi_sq = xi ** 2
        xi_cube = xi ** 3

        # Common denominators
        denom1 = 80.0 * phi ** 2 - 20.0 * phi - 4.0
        denom2 = 32.0 * phi + 8.0
        denom3 = 160.0 * phi ** 2 - 40.0 * phi - 8.0

        # Shape function second derivatives expressions
        shape_functions_second_derivatives_values = np.zeros(6)

        shape_functions_second_derivatives_values[0] = -120.0 * xi_cube / denom3 - \
                                                        48.0 * xi_sq / denom2 + \
                                                        6 * xi * (40.0 * phi + 10.0) / denom3 + \
                                                        2 * (16.0 * phi + 8.0) / denom2

        shape_functions_second_derivatives_values[1] = -12 * L * xi_sq / denom2 + \
                                                        6 * L * xi / denom3 + \
                                                        2 * L / denom2 + \
                                                        20 * xi_cube * (2.0 * L * phi - L) / denom3

        shape_functions_second_derivatives_values[2] = 96.0 * xi_sq / denom2 + \
                                                        2 * (-32.0 * phi - 16.0) / denom2

        shape_functions_second_derivatives_values[3] = 20 * xi_cube * (-4.0 * L * phi - 4.0 * L) / denom3 + \
                                                        6 * xi * (40.0 * L * phi + 8.0 * L) / denom3

        shape_functions_second_derivatives_values[4] = 120.0 * xi_cube / denom3 - \
                                                        48.0 * xi_sq / denom2 + \
                                                        6 * xi * (-40.0 * phi - 10.0) / denom3 + \
                                                        2 * (16.0 * phi + 8.0) / denom2

        shape_functions_second_derivatives_values[5] = 12 * L * xi_sq / denom2 + \
                                                        6 * L * xi / denom3 - \
                                                        2 * L / denom2 + \
                                                        20 * xi_cube * (2.0 * L * phi - L) / denom3
        return shape_functions_second_derivatives_values

    # ------------------------------------------------------------------------------------------------
    def GetThirdDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Precompute powers of xi
        xi_sq = xi ** 2

        # Common denominators
        denom1 = 80.0 * phi ** 2 - 20.0 * phi - 4.0
        denom2 = 32.0 * phi + 8.0
        denom3 = 160.0 * phi ** 2 - 40.0 * phi - 8.0

        # Shape function third derivatives expressions
        shape_functions_third_derivatives_values = np.zeros(6)

        shape_functions_third_derivatives_values[0] = -360.0 * xi_sq / denom3 - \
                                                    96.0 * xi / denom2 + \
                                                    6 * (40.0 * phi + 10.0) / denom3

        shape_functions_third_derivatives_values[1] = -24 * L * xi / denom2 + \
                                                    6 * L / denom3 + \
                                                    60 * xi_sq * (2.0 * L * phi - L) / denom3

        shape_functions_third_derivatives_values[2] = 192.0 * xi / denom2

        shape_functions_third_derivatives_values[3] = 60 * xi_sq * (-4.0 * L * phi - 4.0 * L) / denom3 + \
                                                    6 * (40.0 * L * phi + 8.0 * L) / denom3

        shape_functions_third_derivatives_values[4] = 360.0 * xi_sq / denom3 - \
                                                    96.0 * xi / denom2 + \
                                                    6 * (-40.0 * phi - 10.0) / denom3

        shape_functions_third_derivatives_values[5] = 24 * L * xi / denom2 + \
                                                    6 * L / denom3 + \
                                                    60 * xi_sq * (2.0 * L * phi - L) / denom3
        return shape_functions_third_derivatives_values
    # ------------------------------------------------------------------------------------------------
    def GetN_u(self, xi): # u_0 = N_u * U_e
        # These shape functions are used to interpolate the longitudinal displacements than induce axial forces
        return np.array([0.5 * xi * (xi - 1.0), 
                         (1.0 - xi**2),
                         0.5 * xi * (xi + 1.0)])
    # ------------------------------------------------------------------------------------------------
    def GetN_u_Derivatives(self, xi):
        return 2.0 / self.Length * np.array([xi - 0.5,
                                             -2.0 * xi,
                                             xi + 0.5])