import json
import numpy as np
import os
import matplotlib.pyplot as plt

class StrainCaseHistory:
    def __init__(self, DataDirectory, load_n_steps = 0):
        with open(DataDirectory, "r") as file:
            data = json.load(file)

        self.matrix_strain = np.array(data["Mean_Value_of_GREEN_LAGRANGE_STRAIN_VECTOR"]["Segment_1"])
        self.matrix_stress = np.array(data["Mean_Value_of_CAUCHY_STRESS_VECTOR"]["Segment_1"])

        if load_n_steps > 0:
            # In case we want to load a sample of all strain points
            self.matrix_strain = self.matrix_strain[0:load_n_steps, :]
            self.matrix_stress = self.matrix_stress[0:load_n_steps, :]

        self.n_steps     = self.matrix_strain.shape[0]
        self.strain_size = self.matrix_strain.shape[1]

        self.echo = 0
        if self.echo > 0:
            print("Case: ", DataDirectory, " has been loaded...")


    def IsotropizeData(self, T_strain, T_stress):
        for step in range(self.n_steps):
            self.matrix_strain[step, :] = T_strain @ self.matrix_strain[step, :]
            self.matrix_stress[step, :] = T_stress @ self.matrix_stress[step, :]


    def PrintData(self, component = 0):
        plt.plot(self.matrix_strain[:, component], self.matrix_stress[:, component], color="r")
        plt.grid()
        plt.xlabel("Strain [-]")
        plt.ylabel("Stress [Pa]")
        plt.show()


    def ComputeInternalWorksVector(self):
        accumulated_work = 0.0
        internal_work_vector = np.zeros(self.n_steps)
        delta_strain = np.array([])
        for step in range(0, self.n_steps):
            if step == 0:
                delta_strain = self.matrix_strain[step, :]
            else:
                delta_strain = self.matrix_strain[step, :] - self.matrix_strain[step - 1, :]
            accumulated_work += np.dot(delta_strain, self.matrix_stress[step, :])
            internal_work_vector[step] = accumulated_work
        return internal_work_vector
