import strain_stress_history_case as s_e_history
import constitutive_law_simulated_history

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMApp
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp

import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
import os

'''
    This script does the following operations:
        1 - Read the RVE strain-stress information precomputed (in the same directory inside a folder)
        2 - Computes an approximation of the anisotropic constitutive matrix C_aniso
        3 - Implements an objective function whose "cost" is the difference of internal work dissipated by the
            RVE and a custom isotropic homogenous CL. In this case the CL is a plane stress elastic CL
        4 - Minimizes the value of the objective function by opimizing the young modulus and poisson ratio until convergence
        5 - Returns the C_aniso, C_iso and its Mapping matrix T.

    Author: A. Cornejo, CIMNE, 2024
'''

# We load the cases...
history_cases = []
number_of_steps = 0             # Number of steps to be loaded, 0 means all
dir = "original_elastic_branch" # original rve with custom CL only elastic branches

history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_01_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_02_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_03_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_04_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_05_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_06_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_07_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_08_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_09_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_10_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_11_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_12_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_13_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_14_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_15_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_16_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_17_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_18_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_19_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_20_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_21_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_22_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_23_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_24_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_25_average_Structure.wall_total.json"), number_of_steps))
history_cases.append(s_e_history.StrainCaseHistory(os.path.join(dir, "Case_26_average_Structure.wall_total.json"), number_of_steps))

# Now we compute C_aniso
stress = []
strain = []
for case in history_cases:
    stress.append(case.matrix_stress)
    strain.append(case.matrix_strain)
stress_stack = np.vstack(stress)
strain_stack = np.vstack(strain)
C_aniso = stress_stack.T @ np.linalg.pinv(strain_stack.T)
C_aniso = 0.5*(C_aniso + C_aniso.T)

def ObjetiveFunction(optimizable_property_list):
    '''
    This function reads a set of RVE simulation cases and returns the difference of
    internal works between the reference strain-stress pairs and the modelled by a 
    custom CL

    optimizable_property_list = [
                                YOUNG_MODULUS,
                                POISSON_RATIO
                                ]
    returns a scalar value of the error in terms of internal work.
    '''

    # Now we prepare the CL data...
    current_model = KM.Model()
    model_part = current_model.CreateModelPart("test")
    process_info = model_part.ProcessInfo

    # We extract the nodes from one FE of the macro model geometry
    node1 = model_part.CreateNewNode(1, 0.0,   0.0, 0.0)
    node2 = model_part.CreateNewNode(2, 0.006, 0.0, 0.0)
    node3 = model_part.CreateNewNode(3, 0.006, 0.006, 0.0)
    node4 = model_part.CreateNewNode(4, 0.0,   0.006, 0.0)
    geometry = KM.Quadrilateral2D4(node1, node2, node3, node4)

    cl = SMApp.LinearElasticPlaneStress2DLaw()
    # cl = CLApp.DamageDPlusDMinusPlaneStressMasonry2DLaw()
    cl_pointer = cl.Clone()

    properties = KM.Properties(1)
    properties.SetValue(KM.YOUNG_MODULUS, optimizable_property_list[0] * 1.0e9) # NOTE: we multiply to pass to Pa
    properties.SetValue(KM.POISSON_RATIO, optimizable_property_list[1])
    properties.SetValue(KM.DENSITY, 2400.0)
    properties.SetValue(KM.THICKNESS, 0.311)
    properties.SetValue(KM.CONSTITUTIVE_LAW, cl)

    # We compute the norm of the Wint - Wint_simulation and store it here
    norm_internal_work_difference = np.zeros(len(history_cases))
    for index, case in enumerate(history_cases):

        micro_internal_work_vector = case.ComputeInternalWorksVector()

        # Let's compute the simulated int work...
        strain_history    = case.matrix_strain
        simulated_CL_resp = constitutive_law_simulated_history.ConstitutiveLawReponse(cl_pointer, strain_history, geometry, properties)
        simulated_CL_resp.CalculateStressHistoryResponse()
        simulated_internal_work_vector = simulated_CL_resp.ComputeInternalWorksVector()

        # store the norm of the difference
        internal_work_error = simulated_internal_work_vector - micro_internal_work_vector

        norm_internal_work_difference[index] = internal_work_error[-1]

    # we return the sqrt of the sumation of squared errors of each case
    return np.linalg.norm(norm_internal_work_difference)

# Let's optimize E and nu
initial_guess = [
                3, # E, in GPa
                0.12,  # nu
                ]

my_bounds = [
            (0.1, 10), # E, in GPa
            (0.01, 0.48),  # nu
            ]

result = sp.optimize.minimize(
                        ObjetiveFunction,
                        initial_guess,
                        bounds = my_bounds
                        )
# Print the results
print("\n##############################################################")
print("Optimal value of [E , nu]:", result.x, " [GPa, -]")
print("##############################################################")

print("\n C_aniso: \n", C_aniso)

C_aniso_sqrt = sp.linalg.sqrtm(C_aniso)

# Let's build the isotropic constitutive matrix
E  = result.x[0] * 1.0e9
nu = result.x[1]
d11 = E / (1.0 - nu**2)
d22 = d11
d12 = nu * d11
d33 = E / (2.0 * (1.0 + nu))
C_iso = np.array([
                 [d11, d12, 0.0],
                 [d12, d22, 0],
                 [0.0,0.0,d33]
                 ])

print("\n C_iso: \n", C_iso)
C_iso_sqrt = sp.linalg.sqrtm(C_iso)

T = np.linalg.inv(C_iso_sqrt) @ C_aniso_sqrt
print("\n T: \n", T)

'''
Output:

    ####################################################
    Optimal value of [E , nu]: [4.46701076 0.21639363]  [GPa, -]
    ####################################################

    C_aniso:
    [[5.44213174e+09 8.33854598e+08 4.99364680e+02]
    [8.33854598e+08 4.29133230e+09 6.70074828e+02]
    [4.99364680e+02 6.70074828e+02 1.70685281e+09]]

    C_iso:
    [[4.68645989e+09 1.01412006e+09 0.00000000e+00]
    [1.01412006e+09 4.68645989e+09 0.00000000e+00]
    [0.00000000e+00 0.00000000e+00 1.83616991e+09]]

    T:
    [[ 1.08377292e+00 -1.68525687e-02  5.00845491e-08]
    [-3.03588719e-02  9.60420898e-01  8.35917693e-08]
    [ 9.40748901e-08  1.41461255e-07  9.64143334e-01]]
'''






