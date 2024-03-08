import numpy as np

# Solves K*u = f
def SolveSystem(K, f): # returns u
    return np.dot(np.linalg.inv(K), f)

def ComputeReactions(K, u):
    return np.dot(K, u)

def BuildGlobalSystemContribution(Kglobal, Klocal, ElementDoFList):
    local_size = np.shape(Klocal)[0]
    for row in range(local_size):
        for col in range(local_size):
            Kglobal[ElementDoFList[row], ElementDoFList[col]] += Klocal[row, col]

def BuildSystem(Kglobal, Elements): # Elements is a list of Element
    for element in Elements:
        K_local  = element.CalculateStiffnessMatrix()
        dof_list = element.GetDoFList()
        BuildGlobalSystemContribution(Kglobal, K_local, dof_list)

def InitializeGlobalSystem(NumberOfNodes, dimension): # returns null global K and f 
    global_size = dimension*NumberOfNodes
    return np.zeros((global_size, global_size)), np.zeros(global_size)