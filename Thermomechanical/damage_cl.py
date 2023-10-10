import matplotlib.pyplot as plt
import math
from sympy import *
import scienceplots


"""
This script implements a 1-d damage CL to do some trials when changing the material properties according to Temperature
"""
def ComputeDamageLinear(UniaxialStress, Threshold, A):
    damage = (1.0 - Threshold / UniaxialStress) / (1.0 + A)
    if damage >= 0.9999999:
        damage = 0.9999999
    return damage

def ComputeDamageExponential(UniaxialStress, Threshold, A):
    damage = 1.0 - (Threshold / UniaxialStress) * exp(A * (1.0 - UniaxialStress / Threshold))
    if damage >= 0.9999999:
        damage = 0.9999999
    return damage

def ComputeAParameterLinear(InitialThreshold, Gf, L, Young):
    return -InitialThreshold**2 / (2.0 * Young * Gf / L)

def ComputeAParameterExponential(InitialThreshold, Gf, L, Young):
    return 1.0 / (Gf * Young / (L * InitialThreshold**2) - 0.5)
def ComputeDamage(UniaxialStress, Threshold, Gf, L, Young, Softening):
    if (Softening ==  "Linear"):
        A = ComputeAParameterLinear(Threshold, Gf, L, Young)
        return ComputeDamageLinear(UniaxialStress, Threshold, A)
    else: # Exponential
        A = ComputeAParameterExponential(Threshold, Gf, L, Young)
        return ComputeDamageExponential(UniaxialStress, Threshold, A)

class DamageCL():
    def __init__(self, E, G, fy, l_char):
        self.E = E
        self.G = G
        self.threshold = fy
        self.characteristic_length = l_char
        self.initial_threshold = fy
        self.damage = 0.0


    def CalculateStress(self, strain):
        effective_stress = self.E * strain
        is_damage = False
        if effective_stress >= self.threshold:
            self.threshold = effective_stress
            self.damage = ComputeDamage(effective_stress, self.initial_threshold, self.G, self.characteristic_length, self.E, "Exponential")
            is_damage = True
            return effective_stress * (1.0 - self.damage), is_damage
        else:
            return effective_stress * (1.0 - self.damage), is_damage
# End Class


def CreateStrainHistory(StrainList, max_strain, steps):
    initial_strain = StrainList[0]
    strain_increment = max_strain / steps
    for step in range(steps + 1):
        initial_strain += strain_increment
        # if step > steps/5:
        #     initial_strain -= strain_increment
        # else:
        #     initial_strain += strain_increment
        StrainList.append(initial_strain)


# ===========================================================
E = 30e9      # Pa
G = 250       # J/m2
fy = 2e6      # Pa
l_char = 0.1  # m
strain_history = [0.0] # Adimensional, strain history array

# E = 15e9      # Pa
# G = 250*0.7       # J/m2
# fy = 1e6      # Pa
# l_char = 0.1  # m
# strain_history = [0.0] # Adimensional, strain history array
# ===========================================================

CreateStrainHistory(strain_history, 1.0e-2, 1000)

cl = DamageCL(E, G, fy, l_char)

stress_list = []
predicted_g = 0.5/E*fy**2
G_history = []
damage_history = []

count = 0
# Strain increment process.......
for count, strain in enumerate(strain_history):
    # Here the cl computes stresses at each step
    if count == 200:
        cl.E *= 0.5
        cl.threshold *= 0.5
        cl.initial_threshold *= 0.5
        cl.G *= 0.7
        """
        When modifying the E, we should maintain A parameters to be smooth?
        """
        # cl.G *= 0.5 # the A parameters must be constant??
        # cl.Initial_threshold *= 0.5
        pass
    current_stress, is_damage = cl.CalculateStress(strain)
    stress_list.append(current_stress)
    strain_increment = abs(strain_history[count] - strain_history[count-1])
    if stress_list[count] >= 0.0:
        if is_damage:
            predicted_g += strain_increment*abs(stress_list[count])
        G_history.append(predicted_g*l_char)
        damage_history.append(cl.damage)
# End loading process

print("\n    The simulated G is: ", str(predicted_g*l_char), " J/m2.")
print("\n    The simulated G has a % relative error of: ", str((predicted_g*l_char-G)/100.0), "%.")


# Write a txt file
# Open the file in write mode
with open("data.txt", "w") as file:
    # Use a loop to generate and write numbers to the file
    for i in range(1, len(strain_history)):
        
        file.write("{:.2e}".format(strain_history[i]) + "    " + "{:.2e}".format(stress_list[i]) + "\n")  # Write each number followed by a newline character



fig, ax = plt.subplots(2, 2)
plt.style.use(['science', 'high-vis'])

ax[0, 0].plot(strain_history, stress_list, 'b', linewidth=2)
ax[0, 0].set_xlabel("Strain [-]")
ax[0, 0].set_ylabel("Stress [Pa]")
ax[0, 0].grid()

ax[0, 1].plot(strain_history, G_history, 'r', linewidth=2)
ax[0, 1].set_xlabel("Strain [-]")
ax[0, 1].set_ylabel("Dissipated energy [J/m2]")
ax[0, 1].grid()

ax[1, 0].plot(strain_history, damage_history, 'g', linewidth=2)
ax[1, 0].set_xlabel("Strain [-]")
ax[1, 0].set_ylabel("Damage [-]")
ax[1, 0].grid()

plt.show()



