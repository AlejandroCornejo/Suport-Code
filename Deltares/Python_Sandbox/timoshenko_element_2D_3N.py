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
        self.IntegrationOrder = 4
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

        shape_functions_values[0] = (-40.0 * phi ** 2 - 10.0 * phi) / denom1 * xi + (16.0 * phi + 8.0) / denom2 * xi_sq + (40.0 * phi + 10.0) / denom3 * xi_cube + (-4.0) / denom2 * xi_quad + (-6.0) / denom3 * xi_quint
        shape_functions_values[1] = (-L * phi) / denom1 * xi + L / denom2 * xi_sq + L / denom3 * xi_cube + (-L) / denom2 * xi_quad + (2.0 * L * phi - L) / denom3 * xi_quint
        shape_functions_values[2] = 1.0 + (-32.0 * phi - 16.0) / denom2 * xi_sq + 8.0 / denom2 * xi_quad
        shape_functions_values[3] = (-18.0 * L * phi - 2.0 * L) / denom1 * xi + (40.0 * L * phi + 8.0 * L) / denom3 * xi_cube + (-4.0 * L * phi - 4.0 * L) / denom3 * xi_quint
        shape_functions_values[4] = (40.0 * phi ** 2 + 10.0 * phi) / denom1 * xi + (16.0 * phi + 8.0) / denom2 * xi_sq + (-40.0 * phi - 10.0) / denom3 * xi_cube + (-4.0) / denom2 * xi_quad + 6.0 / denom3 * xi_quint
        shape_functions_values[5] = (-L * phi) / denom1 * xi + (-L) / denom2 * xi_sq + L / denom3 * xi_cube + L / denom2 * xi_quad + (2.0 * L * phi - L) / denom3 * xi_quint
        return shape_functions_values

    # ------------------------------------------------------------------------------------------------
    def GetFirstDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Shape function derivatives expressions
        shape_functions_derivatives_values = np.zeros(6)

        shape_functions_derivatives_values[0] = -30.0*xi**4/(160.0*phi**2 - 40.0*phi - 8.0) - 16.0*xi**3/(32.0*phi + 8.0) + 3*xi**2*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*xi*(16.0*phi + 8.0)/(32.0*phi + 8.0) + (-40.0*phi**2 - 10.0*phi)/(80.0*phi**2 - 20.0*phi - 4.0)
        shape_functions_derivatives_values[1] = -L*phi/(80.0*phi**2 - 20.0*phi - 4.0) - 4*L*xi**3/(32.0*phi + 8.0) + 3*L*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) + 2*L*xi/(32.0*phi + 8.0) + 5*xi**4*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_derivatives_values[2] = 32.0*xi**3/(32.0*phi + 8.0) + 2*xi*(-32.0*phi - 16.0)/(32.0*phi + 8.0)
        shape_functions_derivatives_values[3] = 5*xi**4*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 3*xi**2*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + (-18.0*L*phi - 2.0*L)/(80.0*phi**2 - 20.0*phi - 4.0)
        shape_functions_derivatives_values[4] = 30.0*xi**4/(160.0*phi**2 - 40.0*phi - 8.0) - 16.0*xi**3/(32.0*phi + 8.0) + 3*xi**2*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*xi*(16.0*phi + 8.0)/(32.0*phi + 8.0) + (40.0*phi**2 + 10.0*phi)/(80.0*phi**2 - 20.0*phi - 4.0)
        shape_functions_derivatives_values[5] = -L*phi/(80.0*phi**2 - 20.0*phi - 4.0) + 4*L*xi**3/(32.0*phi + 8.0) + 3*L*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 2*L*xi/(32.0*phi + 8.0) + 5*xi**4*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)

        return shape_functions_derivatives_values / self.J
    # ------------------------------------------------------------------------------------------------
    def GetSecondDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        # Shape function second derivatives expressions
        shape_functions_second_derivatives_values = np.zeros(6)

        shape_functions_second_derivatives_values[0] = -120.0*xi**3/(160.0*phi**2 - 40.0*phi - 8.0) - 48.0*xi**2/(32.0*phi + 8.0) + 6*xi*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*(16.0*phi + 8.0)/(32.0*phi + 8.0)
        shape_functions_second_derivatives_values[1] = -12*L*xi**2/(32.0*phi + 8.0) + 6*L*xi/(160.0*phi**2 - 40.0*phi - 8.0) + 2*L/(32.0*phi + 8.0) + 20*xi**3*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_second_derivatives_values[2] = 96.0*xi**2/(32.0*phi + 8.0) + 2*(-32.0*phi - 16.0)/(32.0*phi + 8.0)
        shape_functions_second_derivatives_values[3] = 20*xi**3*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 6*xi*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_second_derivatives_values[4] = 120.0*xi**3/(160.0*phi**2 - 40.0*phi - 8.0) - 48.0*xi**2/(32.0*phi + 8.0) + 6*xi*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) + 2*(16.0*phi + 8.0)/(32.0*phi + 8.0)
        shape_functions_second_derivatives_values[5] = 12*L*xi**2/(32.0*phi + 8.0) + 6*L*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 2*L/(32.0*phi + 8.0) + 20*xi**3*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)

        return shape_functions_second_derivatives_values / self.J**2

    # ------------------------------------------------------------------------------------------------
    def GetThirdDerivativesShapeFunctionsValues(self, xi):
        # Intermediate variables
        phi = self.Phi
        L = self.Length

        shape_functions_third_derivatives_values = np.zeros(6)
        shape_functions_third_derivatives_values[0] = -360.0*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0*xi/(32.0*phi + 8.0) + 6*(40.0*phi + 10.0)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_third_derivatives_values[1] = -24*L*xi/(32.0*phi + 8.0) + 6*L/(160.0*phi**2 - 40.0*phi - 8.0) + 60*xi**2*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)   
        shape_functions_third_derivatives_values[2] = 192.0*xi/(32.0*phi + 8.0)
        shape_functions_third_derivatives_values[3] = 60*xi**2*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) + 6*(40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)       
        shape_functions_third_derivatives_values[4] = 360.0*xi**2/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0*xi/(32.0*phi + 8.0) + 6*(-40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_third_derivatives_values[5] = 24*L*xi/(32.0*phi + 8.0) + 6*L/(160.0*phi**2 - 40.0*phi - 8.0) + 60*xi**2*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
        return shape_functions_third_derivatives_values / self.J**3
    # ------------------------------------------------------------------------------------------------
    def GetFourthDerivativesShapeFunctionsValues(self, xi):
        phi = self.Phi
        L = self.Length

        shape_functions_fourth_derivatives_values = np.zeros(6)
        shape_functions_fourth_derivatives_values[0] = -720.0*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0/(32.0*phi + 8.0)
        shape_functions_fourth_derivatives_values[1] = -24*L/(32.0*phi + 8.0) + 120*xi*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_fourth_derivatives_values[2] = 192.0/(32.0*phi + 8.0)
        shape_functions_fourth_derivatives_values[3] = 120*xi*(-4.0*L*phi - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0)
        shape_functions_fourth_derivatives_values[4] = 720.0*xi/(160.0*phi**2 - 40.0*phi - 8.0) - 96.0/(32.0*phi + 8.0)
        shape_functions_fourth_derivatives_values[5] = 24*L/(32.0*phi + 8.0) + 120*xi*(2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0)
        return shape_functions_fourth_derivatives_values / self.J**4

    # ------------------------------------------------------------------------------------------------
    def GetN_u(self, xi): # u_0 = N_u * U_e
        # These shape functions are used to interpolate the longitudinal displacements than induce axial forces
        return np.array([0.5 * xi * (xi - 1.0), 
                         (1.0 - xi**2),
                         0.5 * xi * (xi + 1.0)])
    # ------------------------------------------------------------------------------------------------
    def GetN_u_Derivatives(self, xi):
        return 1.0 / self.J * np.array([xi - 0.5,
                                        -2.0 * xi,
                                        xi + 0.5])
    # ------------------------------------------------------------------------------------------------
    def GetN_theta(self, xi): # theta = N_theta * U_e
        N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(xi)
        N_third_deriv = self.GetThirdDerivativesShapeFunctionsValues(xi)
        return N_derivatives + self.Phi * self.Length**2 / 12.0 * N_third_deriv # v' + phi*L**2 / 12 * v'''
    # ------------------------------------------------------------------------------------------------
    def GetN_theta_derivatives(self, xi):
        N_second_derivatives = self.GetSecondDerivativesShapeFunctionsValues(xi)
        N_fourth_derivatives = self.GetFourthDerivativesShapeFunctionsValues(xi)
        return N_second_derivatives + self.Phi * self.Length**2 / 12.0 * N_fourth_derivatives
    # ------------------------------------------------------------------------------------------------
    def GetAxialStrain(self, xi, Ue): # Ue complete 9 components [u0_i, v_i, theta_i]
        Nu_derivatives = self.GetN_u_Derivatives(xi)
        return Nu_derivatives[0] * Ue[0] + Nu_derivatives[3] * Ue[3] + Nu_derivatives[6] * Ue[6]
    # ------------------------------------------------------------------------------------------------
    def GetCurvature(self, xi, Ue): # Ue complete 9 components [u0_i, v_i, theta_i]
        N_theta_derivatives = self.GetN_theta_derivatives(xi)
        return N_theta_derivatives[0] * Ue[1] + N_theta_derivatives[1] * Ue[2] + N_theta_derivatives[2] * Ue[4] + N_theta_derivatives[3] * Ue[5] + N_theta_derivatives[4] * Ue[7] + N_theta_derivatives[5] * Ue[8]
    # ------------------------------------------------------------------------------------------------
    def GetShearStrain(self, xi, U_e): # Ue complete 9 components [u0_i, v_i, theta_i]
        N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(xi)
        N_theta = GetN_theta(xi)
        N_s = N_derivatives - N_theta
        return N_s[0] * Ue[1] + N_s[1] * Ue[2] + N_s[2] * Ue[4] + N_s[3] * Ue[5] + N_s[4] * Ue[7] + N_s[5] * Ue[8]

    def CalculateStiffnessMatrix(self, IntegrationOrder = 2):
        K = np.zeros((9, 9))
        integration_point_w = self.GetIntegrationPoints(IntegrationOrder)
        for IP in range(IntegrationOrder):

            global_size_N = np.zeros(9)
            x_ip = integration_point_w[IP, 0]
            w_ip = integration_point_w[IP, 1]

            # The axial contributions...
            N_u_derivatives = self.GetN_u_Derivatives(x_ip)
            global_size_N[0] = N_u_derivatives[0]
            global_size_N[3] = N_u_derivatives[1]
            global_size_N[6] = N_u_derivatives[2]
            K += np.outer(global_size_N, global_size_N) * self.dN_dEl() * w_ip * self.J # we add the Jacobian of the isoparametric space

            # The bending contributions...
            global_size_N = np.zeros(9)
            N_theta_derivatives = self.GetN_theta_derivatives(x_ip)
            global_size_N[1] = N_theta_derivatives[0]
            global_size_N[2] = N_theta_derivatives[1]
            global_size_N[4] = N_theta_derivatives[2]
            global_size_N[5] = N_theta_derivatives[3]
            global_size_N[7] = N_theta_derivatives[4]
            global_size_N[8] = N_theta_derivatives[5]
            K += np.outer(global_size_N, global_size_N) * self.dM_dkappa() * w_ip * self.J

            # The shear contributions...
            N_theta = self.GetN_theta(x_ip)
            N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(x_ip)
            global_size_N[1] = N_derivatives[0] - N_theta[0]
            global_size_N[2] = N_derivatives[1] - N_theta[1]
            global_size_N[4] = N_derivatives[2] - N_theta[2]
            global_size_N[5] = N_derivatives[3] - N_theta[3]
            global_size_N[7] = N_derivatives[4] - N_theta[4]
            global_size_N[8] = N_derivatives[5] - N_theta[5]
            K += np.outer(global_size_N, global_size_N) * self.dV_dgamma() * w_ip * self.J
        return K