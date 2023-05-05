import matplotlib.pyplot as plt
import math
import numpy

def ComputeDamage(UniaxialStress, Threshold, Gf, L, Young, Softening):
    if (Softening ==  "Linear"):
        A = ComputeAParameterLinear(Threshold, Gf, L, Young)
        return ComputeDamageLinear(UniaxialStress, Threshold, A)
    else: # Exponential
        A = ComputeAParameterExponential(Threshold, Gf, L, Young)
        if A <= 0.0:
            raise ValueError("A negative")
        return ComputeDamageExponential(UniaxialStress, Threshold, A)
def ComputeDamageLinear(UniaxialStress, Threshold, A):
    d = (1.0 - Threshold / UniaxialStress) / (1.0 + A)
    if d >= 0.999:
        return 0.999
    else:
        return d
def ComputeDamageExponential(UniaxialStress, Threshold, A):
    d = 1.0 - (Threshold / UniaxialStress) * math.exp(A * (1.0 - UniaxialStress / Threshold))
    if d >= 0.999:
        return 0.999
    else:
        return d
def ComputeAParameterLinear(Threshold, Gf, L, Young):
    return -Threshold**2 / (2.0 * Young * Gf / L)
def ComputeAParameterExponential(Threshold, Gf, L, Young):
    return 1.0 / (Gf * Young / (L * Threshold**2) - 0.5)

# Parameters
E = 2.1e9      # Pa
l_char = 0.1   # m
fy = 500e6     # Pa
Gf = 15000000      # J/m2
softening = "Linear"

initial_strain = 0.0
final_strain   = fy/E*3.0

strain = numpy.linspace(initial_strain, final_strain, 100)
effective_stress = E*strain

plt.plot(strain, effective_stress, "g", label = "Elastic response", linewidth=2)

integrated_stress = []
for stress_step in effective_stress:
    if stress_step >= fy:
        damage = ComputeDamage(stress_step,fy,Gf,l_char,E,softening)
        integrated_stress.append((1.0-damage)*stress_step)
    else:
        integrated_stress.append(stress_step)

plt.plot(strain, integrated_stress, "r", label = "Damage response", linewidth=2)

plt.xlabel('Strain [-]', fontsize = 12)
plt.ylabel('Stress [Pa]', fontsize = 12)
plt.legend(fontsize = 13)
plt.title('Damage case: ' + "E = " + str(E) + "\n" + "Yield: " + str(fy) + "\n" + "Gf = " + str(Gf) + "\n" + "L_char = " + str(l_char), fontsize = 12)
plt.grid()
plt.show()