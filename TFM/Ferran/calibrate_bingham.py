import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import math

"""
En este script leemos un fichero de datos y calibramos una ley de Bingam:
    tau = tau_0 + mu * \dot{gamma}
"""

###########################################
def ReadDataFromFile(name_of_the_file=""):
    x = []
    y = []
    file = open(name_of_the_file, "r")
    contents = file.readlines()
    for line in contents:
        numbers = line.split()
        if len(numbers) > 0:
            x.append(float(numbers[0]))
            y.append(float(numbers[1]))
    file.close()
    return x, y
###########################################



# bounded data reading of th exp
min_shear = 0.0
max_shear = 100.0

# experimental data file (two columns expected, X=Tau; Y=dot(Gamma))
input_experimental_data = "data_0_sp_sample_16.dat"

tau_exp, gamma_exp = ReadDataFromFile(input_experimental_data)

plt.plot(gamma_exp, tau_exp, "D", markersize=2,  label = input_experimental_data, color="red")

# Now we fit a Bingham CL
def Bingham(gamma_dot, tau_0, mu):
    return tau_0 + mu*gamma_dot

calibrated_parameters = sp.optimize.curve_fit(Bingham, gamma_exp, tau_exp, p0=[60, 0.5])

# print(calibrated_parameters)
# Now we print the calibrated function
tau_0 = calibrated_parameters[0][0]
mu    = calibrated_parameters[0][1]

gamma_dot_bingham = np.linspace(0,120,num=150)
tau_calibrated_bingham = Bingham(gamma_dot_bingham, tau_0, mu)

plt.plot(gamma_dot_bingham, tau_calibrated_bingham, 'k-', label='Fit with ' + r'${\tau_0}$=' + "{0:.4e}".format(tau_0).rjust(11) + " and " + r'${\mu}$=' + "{0:.4e}".format(mu).rjust(11))


plt.xlabel(r'$\dot{\gamma}$ [1/s]', fontsize = 15)
plt.ylabel(r'${\tau}$ [Pa]', fontsize = 15)
plt.ylim(0, 160.0001)
plt.xlim(0, 120.0001)
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
plt.legend()
plt.grid()
plt.show()