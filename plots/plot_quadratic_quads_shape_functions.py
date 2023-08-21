import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
3-----6-----2
|           |
|           |
7     8     5
|           |
|           |
0-----4-----1
 '''
 
#################
shape_function_node = 8
#################

# We define the shape functions
def N(xi, eta, node):
    fx1 = 0.5 * ( xi - 1 ) * xi;
    fx2 = 0.5 * ( xi + 1 ) * xi;
    fx3 = 1 - xi * xi;
    fy1 = 0.5 * ( eta - 1 ) * eta;
    fy2 = 0.5 * ( eta + 1 ) * eta;
    fy3 = 1 - eta * eta;
    if node == 0:
        return fx1*fy1
    elif node == 1:
        return fx2*fy1
    elif node == 2:
        return fx2*fy2
    elif node == 3:
        return fx1*fy2
    elif node == 8:
        return fx3*fy3
    elif node == 5:
        return fx2*fy3
    # ... there are more

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
fig.colorbar(surf, ax=ax)

# Show the plot
plt.show()
