import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMApp
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp
import numpy as np
import matplotlib.pyplot as plt

number_steps = 500
strain_size  = 3
strain_numpy   = np.zeros((number_steps, strain_size))

# We create a custom strain history
factor = 1.0e-6
for row in range(strain_numpy.shape[0]):
    strain_numpy[row, :] = [row * factor, 0.0, 0.0]

strain_history = KM.Matrix(strain_numpy)
stress_history = KM.Matrix(np.zeros((number_steps, strain_size)))

# We prepare the aux model
current_model = KM.Model()
model_part = current_model.CreateModelPart("test")
process_info = model_part.ProcessInfo

node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
node2 = model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
geometry = KM.Triangle2D3(node1, node2, node3)

# cl = SMApp.LinearElasticPlaneStress2DLaw()
cl = CLApp.SmallStrainIsotropicDamagePlaneStressVonMises()

properties = KM.Properties(1)
properties.SetValue(KM.YOUNG_MODULUS, 30e9)
properties.SetValue(KM.POISSON_RATIO, 0.0)
properties.SetValue(KM.DENSITY, 2400.0)
properties.SetValue(KM.THICKNESS, 1.0)
properties.SetValue(KM.YIELD_STRESS, 1e6)
properties.SetValue(KM.FRACTURE_ENERGY, 150)
properties.SetValue(CLApp.SOFTENING_TYPE, 1)
properties.SetValue(KM.CONSTITUTIVE_LAW, cl)

stress_history = CLApp.ComputeCauchyStressHistoryUtility().ComputeStressHistory(cl.Clone(), strain_history, geometry, properties)

# plot it...
x = np.zeros(strain_numpy.shape[0])
y = np.zeros(strain_numpy.shape[0])

for step in range(strain_numpy.shape[0]):
    x[step] = strain_history[step, 0]
    y[step] = stress_history[step, 0]

plt.plot(x, y)
plt.grid()
plt.xlabel("Strain [-]")
plt.ylabel("Stress [Pa]")
plt.show()

