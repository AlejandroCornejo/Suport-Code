
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
        self.hidden_layer_size = 200  # number of hidden neurons

        # Weights matrices wij^(1) and wij^(2)
        self.w1 = np.random.randn(self.input_layer_size,  self.hidden_layer_size)  # 1x4
        self.w2 = np.random.randn(self.hidden_layer_size, self.output_layer_size)  # 4x1
        # self.w1 = np.zeros((self.input_layer_size,  self.hidden_layer_size))  # 1x4
        # self.w2 = np.zeros((self.hidden_layer_size, self.output_layer_size))  # 4x1

        # derivative of the cost function J with respect to the weights
        self.dJ_dw1  = np.zeros((self.input_layer_size,  self.hidden_layer_size))
        self.dJ_dw2  = np.zeros((self.hidden_layer_size, self.output_layer_size))

    def PrintWeights(self):
        print("Weights: \n")
        print("W1: ", self.w1, "\n")
        print("W2: ", self.w2, "\n")

    def PrintGradients(self):
        print("Gradients: \n")
        print("dJ_dW1: ", self.dJ_dw1, "\n")
        print("dJ_dW2: ", self.dJ_dw2, "\n")

    def Sigmoid(self, z):
        return 1.0 / (1.0 + np.exp(-z))

    def Forward(self, x):
        z2 = np.dot(x, self.w1)   # Activity:   number_of_examples x number_of_hidden_layers
        a2 = self.Sigmoid(z2)     # Activation: number_of_examples x number_of_hidden_layers
        z3 = np.dot(a2, self.w2)  # Activity: 
        y_hat = self.Sigmoid(z3)  # prediction
        return y_hat

    def AccumulatedError_J(self, y, y_hat):
        # let's print the accumulated error ==>  AccumErr = Sum(0.5*(y-yhat)^2)
        error_y = y - y_hat
        accumulated_error = 0.0
        rows, cols = error_y.shape
        for i in range(rows):
            for j in range(cols):
                accumulated_error += error_y[i,j]**2
        return 0.5*accumulated_error


# We define the raw data
x = np.transpose(np.array(([np.linspace(-10,10,50)]), dtype=float)) # needs to be an array with 2 dimensions
y = x**2

# We scale the data...
x_max = x.max()
y_max = y.max()
x = x / x_max
y = y / y_max
# plot reference data
pl.plot(x, y, "b +")

# now we plot the first estimation...
neural_network = NeuralNetwork()
y_hat = neural_network.Forward(x)
pl.plot(x, y_hat, color="k", linewidth=5)

current_accumulated_error = neural_network.AccumulatedError_J(y, y_hat)

# now we start training...
perturbation  = 1.0e-12  # for computing the grandients

steepest_descent_factor = 1e-2 # for updating the weights

# Convergence parameters
tolerance = 1.0e-11
iteration = 0
max_iter = 2000
old_accumulated_error = 1.0
relative_error = 1.0

# interval of echo
interval_iter_print = 50
print_counter = 0

while relative_error > tolerance and iteration < max_iter:
    iteration += 1
    # compute gradients of the wij^(1)
    rows, cols = neural_network.w1.shape
    for i in range(rows):
        for j in range(cols):
            neural_network.w1[i,j]    += perturbation
            y_hat_perturbed            = neural_network.Forward(x)
            perturbed_accum_error      = neural_network.AccumulatedError_J(y, y_hat_perturbed)
            neural_network.w1[i,j]    -= perturbation
            neural_network.dJ_dw1[i,j] = (perturbed_accum_error - current_accumulated_error) / perturbation

    # compute gradients of the wij^(2)
    rows, cols = neural_network.w2.shape
    for i in range(rows):
        for j in range(cols):
            neural_network.w2[i,j]    += perturbation
            y_hat_perturbed            = neural_network.Forward(x)
            perturbed_accum_error      = neural_network.AccumulatedError_J(y, y_hat_perturbed)
            neural_network.w2[i,j]    -= perturbation
            neural_network.dJ_dw2[i,j] = (perturbed_accum_error - current_accumulated_error) / perturbation

    # now we update the weights ==> Learning...
    neural_network.w1 -= neural_network.dJ_dw1 * steepest_descent_factor
    neural_network.w2 -= neural_network.dJ_dw2 * steepest_descent_factor

    y_hat = neural_network.Forward(x)
    current_accumulated_error = neural_network.AccumulatedError_J(y, y_hat)
    relative_error = np.abs((old_accumulated_error - current_accumulated_error) / current_accumulated_error)
    old_accumulated_error = current_accumulated_error

    print_counter += 1
    if print_counter > interval_iter_print:
        print(" ## Iteration: ", iteration)
        print("    The current relative_error is: ", relative_error)
        print("    The current_accumulated_error is: ", current_accumulated_error, "\n")
        print_counter = 0
        pl.plot(x, y_hat, color="r")

pl.plot(x, y_hat, color="g", linewidth=5)
pl.show()
