import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMApp
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp
import numpy as np
import matplotlib.pyplot as plt


class ConstitutiveLawReponse():
    def __init__(self, ConstitutiveLawPointer, StrainHistory, Geometry, Properties):
        self.pCL = ConstitutiveLawPointer
        self.matrix_strain = StrainHistory
        self.geometry = Geometry
        self.properties = Properties
        self.n_steps = StrainHistory.shape[0]
        self.strain_size = StrainHistory.shape[1]
        self.matrix_stress = np.zeros((self.n_steps, self.strain_size))

    # This method returns a matrix of n_steps and response stress vectors
    def CalculateStressHistoryResponse(self):
        self.matrix_stress = np.array(CLApp.ComputeCauchyStressHistoryUtility().ComputeStressHistory(self.pCL, self.matrix_strain, self.geometry, self.properties))
        return self.matrix_stress

    # This method computes the accumulated internal work
    # NOTE: the stress response must be computed before hand!
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

    # We do the transofrmation from KM.Matrix to np.array to print...
    def PrintData(self, component = 0):
        stress_to_print = np.zeros(self.n_steps)
        for step in range(self.n_steps):
            stress_to_print[step] = self.matrix_stress[step, component]
        plt.plot(self.matrix_strain[:, component], stress_to_print)
        plt.grid()
        plt.xlabel("Strain [-]")
        plt.ylabel("Stress [Pa]")
        plt.show()