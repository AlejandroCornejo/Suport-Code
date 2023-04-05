
'''
https://www.youtube.com/watch?v=bxe2T-V8XRs&list=PLiaHhY2iBX9hdHaRr6b7XevZtgZRa1PoU&index=1
https://www.youtube.com/watch?v=UJwK6jAStmg&list=PLiaHhY2iBX9hdHaRr6b7XevZtgZRa1PoU&index=2


'''

import numpy as np
import matplotlib.pyplot as pl

class NeuralNetwork(object):
    def __init__(self):
        # We define Hyperparameters
        self.input_layer_size  = 1  # one dimensional input
        self.output_layer_size = 1  # one dimensional output
        self.hidden_layer_size = 4  # number of hidden neurons

        # Weights matrices
        self.w1 = np.random.randn(self.input_layer_size,  self.hidden_layer_size)
        self.w2 = np.random.randn(self.hidden_layer_size, self.output_layer_size)

    def Sigmoid(self, z):
        return 1.0 / (1.0+np.exp(-z))

    def Forward(self, x):
        z2 = np.dot(x, self.w1)   # Activity:   number_of_examples x number_of_hidden_layers
        a2 = self.Sigmoid(z2)     # Activation: number_of_examples x number_of_hidden_layers
        z3 = np.dot(a2, self.w2)  # Activity: 
        y_hat = self.Sigmoid(z3)  # prediction
        return y_hat

# We define the raw data
x = np.transpose(np.array(([np.linspace(-10,10,20)]), dtype=float)) # needs to be an array with 2 dimensions
y = x**2

pl.plot(x, y)

# now we plot the estimation...
neural_network = NeuralNetwork()
y_hat = neural_network.Forward(x)

pl.plot(x, y_hat)

# let's print the accumulated error ==>  AccumErr = Sum(0.5*(y-yhat)^2)
error_y = y - y_hat
accumulated_error = 0.0
for component in error_y:
    accumulated_error += component**2
accumulated_error = 0.5*accumulated_error

print("The accumulated error is: ", str(accumulated_error))
pl.show()
