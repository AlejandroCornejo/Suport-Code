

class Node2D():
    def __init__(self, Id, x, y = 0.0):
        self.Id = Id
        self.x = x
        self.y = y
        self.u        = 0.0 # axial displacement
        self.v        = 0.0 # Vertical deflection
        self.rotation = 0.0 # Actual rotation Theta

    def GetDoFList(self): # DoF start with 0
        return (self.Id - 1) * 2, (self.Id - 1) * 2 + 1