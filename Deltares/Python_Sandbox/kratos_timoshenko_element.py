import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

current_model = KM.Model()
model_part = current_model.CreateModelPart("test")
process_info = model_part.ProcessInfo

model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT_X)
model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT_Y)
model_part.AddNodalSolutionStepVariable(KM.ROTATION_Z)

L = 5.0
P = -1e3

node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
node2 = model_part.CreateNewNode(2, L,   0.0, 0.0)
# node2 = model_part.CreateNewNode(2, 0.0,   L, 0.0)

node1.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
node1.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
node1.AddDof(KM.ROTATION_Z,     KM.MOMENT_Z)

node2.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
node2.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
node2.AddDof(KM.ROTATION_Z,     KM.MOMENT_Z)

node2.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.001)

cl = SMA.TimoshenkoBeamElasticConstitutiveLaw()

properties = KM.Properties(1)

properties.SetValue(KM.YOUNG_MODULUS,2.1e9)
properties.SetValue(KM.POISSON_RATIO,0.3)
properties.SetValue(SMA.CROSS_AREA, 116.0e-4)
properties.SetValue(SMA.IZ, 48200.0e-8)
properties.SetValue(SMA.SHEAR_CORRECTION_XY, 5.0 / 6.0)
properties.SetValue(KM.DENSITY, 2400.0)
properties.SetValue(KM.CONSTITUTIVE_LAW, cl)

geometry = KM.Line2D2(node1, node2)


cl_options = KM.Flags()
cl_options.Set(KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
cl_options.Set(KM.ConstitutiveLaw.COMPUTE_STRESS, True)
strain_vector = KM.Vector(3)
stress_vector = KM.Vector(3)
const_matrix = KM.Matrix(3,3)

cl_values = KM.ConstitutiveLawParameters()
cl_values.SetOptions(cl_options)
cl_values.SetStrainVector(strain_vector)
cl_values.SetStressVector(stress_vector)
cl_values.SetConstitutiveMatrix(const_matrix)
cl_values.SetElementGeometry(geometry)
cl_values.SetMaterialProperties(properties)

timoshenko_element = model_part.CreateNewElement("TimoshenkoBeamElement2D2N", 1, geometry, properties)
timoshenko_element.Initialize(process_info)

cl.Check(properties, geometry, process_info)

LHS = KM.Matrix(6,6)
RHS = KM.Vector(6)

timoshenko_element.CalculateLocalSystem(LHS, RHS, process_info)
# print(LHS)

print(RHS)