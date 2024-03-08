
import numpy as np
import matplotlib.pyplot as pl
import math

"""  

Felippa and Oñate, "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability", Archives of Comp. Methods in Eng. (2021) 28:2021-2080

In this file we derive the main equations of 2-noded Timoshenko FE, firstly using a closed stiffness form (K) and then by integration (int BtDB dV)

"""

class TimoshenkoElement2D2N():
    # ------------------------------------------------------------------------------------------------
    def __init__(self, Node1, Node2, E, I, nu, A):
        self.Nodes = [Node1, Node2]
        self.E = E
        self.A = A
        self.As = A * 5.0 / 6.0 # NOTE depends on the cross section
        self.I = I
        self.nu = nu
        self.Length = math.sqrt((Node1.x - Node2.x)**2 + (Node1.y - Node2.y)**2)
        self.G = E / (2.0 * (1.0 + nu))
        self.Phi = 12.0 * E * I / (self.G * self.As * self.Length**2)
        self.IntegrationOrder = 2
        self.J = self.Length * 0.5
    # ------------------------------------------------------------------------------------------------
    def GetDoFList(self):
        dofs = np.zeros(4)
        dofs[0], dofs[1] = self.Nodes[0].GetDoFList()
        dofs[2], dofs[3] = self.Nodes[1].GetDoFList()
        return dofs
    # ------------------------------------------------------------------------------------------------
    def dN_dEl(self): # This method return the derivative of the axial force with respect to the axial strain El
        return self.E * self.A
    # ------------------------------------------------------------------------------------------------
    def dM_dkappa(self): # This method return the derivative of the Moment with respect to the cuvature kappa
        return self.E * self.I
    # ------------------------------------------------------------------------------------------------
    def dV_dgamma(self): # This method return the derivative of the shear V with respect to the shear strain gamma
        return self.G * self.As
    # ------------------------------------------------------------------------------------------------
    def GetIntegrationPoints(self, IntegrationOrder = 3):
        if IntegrationOrder == 1:
            return np.matrix([0, 2.0]) # Coord, weight
        elif IntegrationOrder == 2: # enough to integrate exact
            return np.matrix([[-1.0/math.sqrt(3.0), 1.0],
                               [1.0/math.sqrt(3.0), 1.0]])
        elif IntegrationOrder == 3:
            return np.matrix([[-math.sqrt(3.0/5.0), 5.0 / 9.0],
                               [0.0               , 8.0 / 9.0],
                               [math.sqrt(3.0/5.0), 5.0 / 9.0]])

    # ------------------------------------------------------------------------------------------------
    def CalculateDirectK(self):
        L = self.Length
        factor = self.E * self.I / (L**3 * (1.0 + self.Phi))
        K = np.matrix([[12.0,            6.0*L,                  -12.0,              6.0*L],
                       [6*L,   L**2*(4.0 + self.Phi),  -6.0*L,   L**2*(2.0-self.Phi)],
                       [-12.0,           -6.0*L,                 12.0,               -6.0*L],
                       [6*L,   L**2*(2.0-self.Phi),    -6.0*L,   L**2*(4.0+self.Phi)]])
        return factor*K
    # ------------------------------------------------------------------------------------------------

    def GetShapeFunctionsValues(self, xi): # xi is the natural coordinate
        # These shape functions are used to interpolate the y deflection "v"
        shape_functions_values = np.zeros(4)
        one_plus_phi = 1.0 + self.Phi
        xi_square = xi**2
        shape_functions_values[0] = (xi - 1.0) * (xi + xi_square - 2.0 * one_plus_phi) / (4.0 * one_plus_phi) # sign was wrong
        shape_functions_values[1] = (1.0 - xi_square) * (1.0 - xi + self.Phi) * self.Length / (8.0 * one_plus_phi)
        shape_functions_values[2] = (1.0 + xi) * (xi - xi_square + 2.0 * one_plus_phi) / (4.0 * one_plus_phi)
        shape_functions_values[3] = (xi_square - 1.0) * (1.0 + xi + self.Phi) * self.Length / (8.0 * one_plus_phi) # sign was wrong
        return shape_functions_values
    # ------------------------------------------------------------------------------------------------
    def GetFirstDerivativesShapeFunctionsValues(self, xi): # xi is the natural coordinate
        shape_functions_derivatives_values = np.zeros(4)
        one_plus_phi = 1.0 + self.Phi
        xi_square = xi**2
        shape_functions_derivatives_values[0] = (-6.0 + 6.0 * xi_square - 4.0 * self.Phi) / (4.0 * one_plus_phi*self.Length)
        shape_functions_derivatives_values[1] = (-1.0 + 3.0 * xi_square - 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi)
        shape_functions_derivatives_values[2] = (6.0 - 6.0 * xi_square + 4.0 * self.Phi) / (4.0 * one_plus_phi * self.Length)
        shape_functions_derivatives_values[3] = (-1.0 + 3.0 * xi_square + 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi)
        return shape_functions_derivatives_values
    # ------------------------------------------------------------------------------------------------
    def GetSecondDerivativesShapeFunctionsValues(self, xi): # xi is the natural coordinate
        shape_functions_second_derivatives_values = np.zeros(4)
        one_plus_phi = 1.0 + self.Phi
        xi_square = xi**2
        shape_functions_second_derivatives_values[0] = 6.0 * xi / (one_plus_phi * self.Length**2.0)
        shape_functions_second_derivatives_values[1] = (-1.0 + 3.0 * xi - self.Phi) / (one_plus_phi * self.Length)
        shape_functions_second_derivatives_values[2] = -6.0 * xi / (one_plus_phi * self.Length**2)
        shape_functions_second_derivatives_values[3] = (1.0 + 3.0 * xi + self.Phi) / (one_plus_phi * self.Length)
        return shape_functions_second_derivatives_values
    # ------------------------------------------------------------------------------------------------
    def GetThirdDerivativesShapeFunctionsValues(self, xi): # xi is the natural coordinate
        shape_functions_third_derivatives_values = np.zeros(4)
        one_plus_phi = 1.0 + self.Phi
        xi_square = xi**2
        shape_functions_third_derivatives_values[0] = 12.0  / (one_plus_phi * self.Length**3)
        shape_functions_third_derivatives_values[1] = 6.0   / (one_plus_phi * self.Length**2)
        shape_functions_third_derivatives_values[2] = -12.0 / (one_plus_phi * self.Length**3)
        shape_functions_third_derivatives_values[3] = 6.0   / (one_plus_phi * self.Length**2)
        return shape_functions_third_derivatives_values
    # ------------------------------------------------------------------------------------------------
    def GetCurvature(self, xi, U_e): # kappa = N_theta_derivative * U_e
        N_theta_derivatives = self.GetN_theta_derivatives(xi)
        return np.dot(N_theta_derivatives, U_e)

    # ------------------------------------------------------------------------------------------------
    def GetShearStrain(self, xi, U_e): # kappa = N_theta_derivative * U_e
        N_theta = self.GetN_theta(xi)
        N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(xi)
        return np.dot(N_derivatives - N_theta, U_e)
    # ------------------------------------------------------------------------------------------------
    def GetN_theta(self, xi): # theta = N_theta * U_e
        # These shape functions are used to interpolate the total rotation Theta
        one_plus_phi = 1.0 + self.Phi
        return np.array([(3.0 * xi**2 - 3.0) / (2.0 * one_plus_phi * self.Length), 
                         (xi - 1.0)*(1.0 + 3.0 * xi - 2.0 * self.Phi) / (4.0 * one_plus_phi),
                         (3.0 - 3.0 * xi**2.0) / (2 * one_plus_phi * self.Length),
                         (1.0 + xi)*(3.0 * xi - 1.0 + 2.0 * self.Phi) / (4.0 * one_plus_phi)])
    # ------------------------------------------------------------------------------------------------
    def GetN_theta_derivatives(self, xi):
        one_plus_phi = 1.0 + self.Phi
        return (2.0 / self.Length) * np.array([3.0 * xi / (one_plus_phi * self.Length), 
                         (-0.5 * self.Phi + 1.5 * xi - 0.5) / (self.Phi + 1.0),
                         (-3.0 * xi) / (self.Length * self.Phi + self.Length),
                         (0.5*self.Phi + 1.5 * xi + 0.5) / (self.Phi + 1.0)])
    # ------------------------------------------------------------------------------------------------
    def GetN_u(self, xi): # u_0 = N_u * U_e
        # These shape functions are used to interpolate the longitudinal displacements than induce axial forces
        return np.array([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)])
    # ------------------------------------------------------------------------------------------------
    def GetN_u_Derivatives(self, xi):
        return 2.0 / self.Length * np.array([-0.5, 0.5])
    # ------------------------------------------------------------------------------------------------

    def CalculateStiffnessMatrix(self, IntegrationOrder = 2):
        K = np.zeros((6, 6))

        integration_point_w = self.GetIntegrationPoints(IntegrationOrder)
        # Integration loop
        for IP in range(IntegrationOrder):
            global_size_N = np.zeros(6)
            x_ip = integration_point_w[IP, 0]
            w_ip = integration_point_w[IP, 1]

            # The axial contributions...
            N_u_derivatives = self.GetN_u_Derivatives(x_ip)
            global_size_N[0] = N_u_derivatives[0]
            global_size_N[3] = N_u_derivatives[1]
            K += np.outer(global_size_N, global_size_N) * self.dN_dEl() * w_ip * self.J # we add the Jacobian of the isoparametric space

            # The bending contributions...
            global_size_N = np.zeros(6)
            N_theta_derivatives = self.GetN_theta_derivatives(x_ip)
            global_size_N[1] = N_theta_derivatives[0]
            global_size_N[2] = N_theta_derivatives[1]
            global_size_N[4] = N_theta_derivatives[2]
            global_size_N[5] = N_theta_derivatives[3]
            K += np.outer(global_size_N, global_size_N) * self.dM_dkappa() * w_ip * self.J

            # The shear contributions...
            global_size_N = np.zeros(6)
            N_theta = self.GetN_theta(x_ip)
            N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(x_ip)
            global_size_N[1] = N_derivatives[0] - N_theta[0]
            global_size_N[2] = N_derivatives[1] - N_theta[1]
            global_size_N[4] = N_derivatives[2] - N_theta[2]
            global_size_N[5] = N_derivatives[3] - N_theta[3]
            K += np.outer(global_size_N, global_size_N) * self.dV_dgamma() * w_ip * self.J
        return K
    # ------------------------------------------------------------------------------------------------
    def CalculateInternalForcesVector(self, IntegrationOrder = 2)
    # ------------------------------------------------------------------------------------------------

    def ApplyBoundaryConditionsToK(self, K, FixedDoFArray):
        # We put all null except the diagonal (block builder)
        size = np.shape(K)[0]
        K_bc = np.copy(K)
        for dof in FixedDoFArray:
            K_bc[dof, :]   = np.zeros(size)
            K_bc[:, dof]   = np.zeros(size)
            K_bc[dof, dof] = 1.0
        return K_bc

    # PRINT METHODS #
    def PrintStrainKinematics(self, U_e):
        xi = np.linspace(-1, 1, 500)
        kappa = np.zeros(xi.size)
        gamma = np.zeros(xi.size)
        counter = 0
        for x in xi:
            kappa[counter] = self.GetCurvature(x, U_e)
            gamma[counter] = self.GetShearStrain(x, U_e)
            counter += 1
        pl.plot(xi, kappa, label="Curvature")
        pl.plot(xi, gamma, label="Shear Strain")
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintShapeFunctions(self):
        xi = np.linspace(-1, 1, 500)
        N1 = np.zeros(xi.size)
        N2 = np.zeros(xi.size)
        N3 = np.zeros(xi.size)
        N4 = np.zeros(xi.size)
        counter = 0
        for x in xi:
            shape_functions_values = self.GetShapeFunctionsValues(x)
            N1[counter] = shape_functions_values[0]
            N2[counter] = shape_functions_values[1]
            N3[counter] = shape_functions_values[2]
            N4[counter] = shape_functions_values[3]
            counter += 1
        pl.plot(xi, N1, label="N1")
        pl.plot(xi, N2, label="N2")
        pl.plot(xi, N3, label="N3")
        pl.plot(xi, N4, label="N4")
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintDeflectionCurveFromNodalValues(self, U_e): # U_e = {v1, rot1, v2, rot2}
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1, 1, 500)
        deflection = np.zeros(xi.size)
        counter = 0
        for x in xi:
            shape_functions_values = self.GetShapeFunctionsValues(x)
            deflection[counter] = np.dot(shape_functions_values, U_e) # v = N·u
            counter += 1
        pl.plot(xi * self.Length / 2 + x_c, deflection, label="v(x)", linewidth=2) # transformed to X
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintRotationCurveFromNodalValues(self, U_e): # U_e = {v1, rot1, v2, rot2}
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1, 1, 500)
        rotation = np.zeros(xi.size)
        counter = 0
        one_plus_phi = 1.0 + self.Phi
        for x in xi:
            N_theta = self.GetN_theta(x)
            rotation[counter] = np.dot(N_theta, U_e)
            counter += 1
        pl.plot(xi * self.Length / 2. + x_c, rotation, label="Rotation(x) / N_theta") # transformed to X
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintBendingMomentsFromNodalValues(self, Ue):
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1., 1., 500)
        M = np.zeros(xi.size)
        counter = 0
        for x in xi:
            kappa = self.GetCurvature(x, Ue)
            M[counter] = -self.E * self.I * kappa # the spanish convention is positive when tension in lower fibre
            counter += 1
        pl.plot(xi * self.Length / 2. + x_c, M, label="M(x)") # transformed to X
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintShearForceFromNodalValues(self, Ue):
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1.0, 1.0, 500)
        V = np.zeros(xi.size)
        counter = 0
        for x in xi:
            gamma_xy = self.GetShearStrain(x, Ue)
            V[counter] = self.G * self.As * gamma_xy
            counter += 1
        pl.plot(xi * self.Length / 2. + x_c, V, label="V(x)") # transformed to X
        pl.grid()
        # y_ticks = np.arange(-1100, 0, 100)
        # pl.yticks(y_ticks, fontsize = 12)
        pl.legend()
        pl.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        pl.show()