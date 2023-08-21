import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
             v
             ^
             |
       4-----------3
       |     |     |
       |     |     |
       |     +---- | --> u
       |           |
       |           |
       1-----------2
 '''
 
#################
shape_function_node = 3
#################

# We define the shape functions
def N(xi, eta, node):
    if node == 1:
        return 0.25*(1.0-xi)*(1.0-eta)
    elif node == 2:
        return 0.25*(1.0+xi)*(1.0-eta)
    elif node == 3:
        return 0.25*(1.0+xi)*(1.0+eta)
    elif node == 4:
        return 0.25*(1.0-xi)*(1.0+eta)

# Generate x and y values
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(x, y)

# Calculate the function values for each (x, y) pair
Z = N(X, Y, shape_function_node)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
surf = ax.plot_surface(X, Y, Z, cmap='viridis')

# Add labels and title
ax.set_xlabel('xi')
ax.set_ylabel('eta')
ax.set_zlabel('N(xi.eta)')
ax.set_title('Shape function of node ' + str(shape_function_node), fontsize=16)

# Add colorbar

fig.colorbar(surf, ax=ax, ticks=np.linspace(0.1, 1.1, 11))

# Show the plot
plt.show()
