
import numpy as np
import matplotlib.pyplot as pl
import math

"""  

Felippa and Oñate, "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability", Archives of Comp. Methods in Eng. (2021) 28:2021-2080

In this file we derive the main equations of 2-noded Timoshenko FE, firstly using a closed stiffness form (K) and then by integration (int BtDB dV)

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
- N_theta   : Rotation shape functions
- N         : Transverse displacement shape functions

"""

class TimoshenkoElement2D2N():
    # ------------------------------------------------------------------------------------------------
    def __init__(self, Node1, Node2, E, I, nu, A):
        self.Nodes = [Node1, Node2]
        self.E      = E
        self.A      = A
        self.As     = A * 5.0 / 6.0 # NOTE depends on the cross section
        self.I      = I
        self.nu     = nu
        self.Length = math.sqrt((Node1.x - Node2.x)**2 + (Node1.y - Node2.y)**2)
        self.alpha  = self.GetAlphaAngle()
        self.G      = E / (2.0 * (1.0 + nu))
        self.Phi    = 12.0 * E * I / (self.G * self.As * self.Length**2)
        self.IntegrationOrder = 2
        self.J      = self.Length * 0.5
    # ------------------------------------------------------------------------------------------------
    def GetNumberOfNodes(self):
        return len(self.Nodes)
    # ------------------------------------------------------------------------------------------------
    def GetDoFList(self):
        dofs = np.zeros(4)
        for i in range(self.GetNumberOfNodes()):
            dofs[2*i], dofs[2*i+1] = self.Nodes[i].GetDoFList()
        return dofs
    # ------------------------------------------------------------------------------------------------
    def GetAlphaAngle(self):
        node_1 = self.Nodes[0]
        node_2 = []
        if (self.GetNumberOfNodes() == 2):
            node_2 = self.Nodes[1]
        else:
            node_2 = self.Nodes[2]
        delta_x = node_2.x - node_1.x
        delta_y = node_2.y - node_1.y
        if abs(delta_x) > 0.0:
            return math.atan((node_2.y - node_1.y) / (node_2.x - node_1.x))
        else:
            if delta_y > 0.0:
                return 0.5 * math.pi
            else:
                return -0.5 * math.pi
    # ------------------------------------------------------------------------------------------------
    def CalculateRotationMatrix(self):
        T = np.zeros((3, 3))
        c = math.cos(self.alpha)
        s = math.sin(self.alpha)
        T[0,0] = c
        T[0,1] = -s
        T[1,0] = s
        T[1,1] = c
        T[2,2] = 1.0
        return T
    # ------------------------------------------------------------------------------------------------
    def RotateK(self, K): # we rotate and return, from local to global
        # Kaa = K[:3, :3]  Kab = K[:3, 3:]   Kba = K[3:, :3]   Kbb = K[3:, 3:]
        T = self.CalculateRotationMatrix()
        Tt = np.transpose(T)
        K[:3, :3] = np.dot(T, np.dot(K[:3, :3], Tt))
        K[:3, 3:] = np.dot(T, np.dot(K[:3, 3:], Tt))
        K[3:, :3] = np.dot(T, np.dot(K[3:, :3], Tt))
        K[3:, 3:] = np.dot(T, np.dot(K[3:, 3:], Tt))
        return K
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
        elif IntegrationOrder == 4:
            return np.matrix([[ -0.861136, 0.347855],
                               [-0.339981,0.652145],
                               [ 0.339981,0.652145],
                               [ 0.861136, 0.347855]])

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
    def GetAxialStrain(self, xi, Ue): # Ue complete 6 components
        Nu_derivatives = self.GetN_u_Derivatives(xi)
        return Nu_derivatives[0] * Ue[0] + Nu_derivatives[1] * Ue[3]
    # ------------------------------------------------------------------------------------------------
    def GetCurvature(self, xi, U_e): # kappa = N_theta_derivative * U_e
        N_theta_derivatives = self.GetN_theta_derivatives(xi)
        # np.dot(N_theta_derivatives, U_e) in global size
        return N_theta_derivatives[0] * U_e[1] + N_theta_derivatives[1] * U_e[2] + N_theta_derivatives[2] * U_e[4] + N_theta_derivatives[3] * U_e[5]

    # ------------------------------------------------------------------------------------------------
    def GetShearStrain(self, xi, U_e): # kappa = N_theta_derivative * U_e, Ue only v and theta
        N_theta = self.GetN_theta(xi)
        N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(xi)
        N_s = N_derivatives - N_theta
        # np.dot(N_derivatives - N_theta, U_e), U_e) in global size
        return N_s[0] * U_e[1] + N_s[1] * U_e[2] + N_s[2] * U_e[4] + N_s[3] * U_e[5]
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
    def CalculateInternalForcesVector(self, Ue, IntegrationOrder = 2):
        F_int = np.zeros(6)

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
            N = self.GetAxialStrain(x_ip, Ue) * self.dN_dEl()
            F_int += global_size_N * N * w_ip * self.J # we add the Jacobian of the isoparametric space

            # The bending contributions...
            global_size_N = np.zeros(6)
            N_theta_derivatives = self.GetN_theta_derivatives(x_ip)
            global_size_N[1] = N_theta_derivatives[0]
            global_size_N[2] = N_theta_derivatives[1]
            global_size_N[4] = N_theta_derivatives[2]
            global_size_N[5] = N_theta_derivatives[3]
            M = self.GetCurvature(x_ip, Ue) * self.dM_dkappa()
            F_int += global_size_N * M * w_ip * self.J

            # The shear contributions...
            global_size_N = np.zeros(6)
            N_theta = self.GetN_theta(x_ip)
            N_derivatives = self.GetFirstDerivativesShapeFunctionsValues(x_ip)
            global_size_N[1] = N_derivatives[0] - N_theta[0]
            global_size_N[2] = N_derivatives[1] - N_theta[1]
            global_size_N[4] = N_derivatives[2] - N_theta[2]
            global_size_N[5] = N_derivatives[3] - N_theta[3]
            V = self.GetShearStrain(x_ip, Ue) * self.dV_dgamma()
            F_int += global_size_N * V * w_ip * self.J
        return F_int
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
    def PrintStrainKinematics(self, U_e): # U_e global size (6)
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
    def PrintDeflectionCurveFromNodalValues(self, U_e): # U_e = {u01, v1, rot1, u02, v2, rot2}
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1, 1, 500)
        deflection = np.zeros(xi.size)
        counter = 0
        global_size_shape_functions = np.zeros(6)
        for x in xi:
            shape_functions_values = self.GetShapeFunctionsValues(x)
            global_size_shape_functions[1] = shape_functions_values[0]
            global_size_shape_functions[2] = shape_functions_values[1]
            global_size_shape_functions[4] = shape_functions_values[2]
            global_size_shape_functions[5] = shape_functions_values[3]
            deflection[counter] = np.dot(global_size_shape_functions, U_e) # v = N·u
            counter += 1
        pl.plot(xi * self.Length / 2 + x_c, deflection, label="v(x)", linewidth=2) # transformed to X
        pl.grid()
        pl.legend()
        pl.show()
    # ------------------------------------------------------------------------------------------------
    def PrintRotationCurveFromNodalValues(self, U_e): # U_e = {u01, v1, rot1, u02, v2, rot2}
        x_c = 0.5*(self.Nodes[0].x + self.Nodes[1].x)
        xi = np.linspace(-1, 1, 500)
        rotation = np.zeros(xi.size)
        # rotation_trial = np.zeros(xi.size)
        counter = 0
        one_plus_phi = 1.0 + self.Phi
        global_size_shape_functions = np.zeros(6)
        for x in xi:
            N_theta = self.GetN_theta(x)
            global_size_shape_functions[1] = N_theta[0]
            global_size_shape_functions[2] = N_theta[1]
            global_size_shape_functions[4] = N_theta[2]
            global_size_shape_functions[5] = N_theta[3]
            rotation[counter] = np.dot(global_size_shape_functions, U_e)

            # N_deriv = self.GetFirstDerivativesShapeFunctionsValues(x)
            # N_third_deriv = self.GetThirdDerivativesShapeFunctionsValues(x)
            # global_size_shape_functions[1] = N_deriv[0]
            # global_size_shape_functions[2] = N_deriv[1]
            # global_size_shape_functions[4] = N_deriv[2]
            # global_size_shape_functions[5] = N_deriv[3]
            # rotation_trial[counter] = np.dot(global_size_shape_functions, U_e)
            # global_size_shape_functions[1] = N_third_deriv[0]
            # global_size_shape_functions[2] = N_third_deriv[1]
            # global_size_shape_functions[4] = N_third_deriv[2]
            # global_size_shape_functions[5] = N_third_deriv[3]
            # rotation_trial[counter] += np.dot(global_size_shape_functions, U_e) * self.Phi * self.Length**2 / 12.0

            counter += 1
        pl.plot(xi * self.Length / 2. + x_c, rotation, label="Rotation(x) / N_theta") # transformed to X
        # pl.plot(xi * self.Length / 2. + x_c, rotation_trial, label="Rotation (x) / Using DEq") # transformed to X

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