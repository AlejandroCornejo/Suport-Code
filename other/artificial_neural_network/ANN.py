
'''
https://www.youtube.com/watch?v=bxe2T-V8XRs&list=PLiaHhY2iBX9hdHaRr6b7XevZtgZRa1PoU&index=1
https://www.youtube.com/watch?v=UJwK6jAStmg&list=PLiaHhY2iBX9hdHaRr6b7XevZtgZRa1PoU&index=2

Artificial neural network with only one hidden layer of "n" neurons.
'''

import numpy as np
import matplotlib.pyplot as pl

class NeuralNetwork(object):
    def __init__(self):
        # We define Hyperparameters
        self.input_layer_size  = 1  # one dimensional input
        self.output_layer_size = 1  # one dimensional output
        self.hidden_layer_size = 25  # number of hidden neurons

        # Weights matrices wij^(1) and wij^(2)
        self.w1 = np.random.randn(self.input_layer_size,  self.hidden_layer_size)  # 1x4
        self.w2 = np.random.randn(self.hidden_layer_size, self.output_layer_size)  # 4x1

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
    
    def Compute_dJdw_finite_differences(self, perturbation, current_accumulated_error):
        rows, cols = self.w1.shape
        for i in range(rows):
            for j in range(cols):
                self.w1[i,j]         += perturbation
                y_hat_perturbed       = self.Forward(x)
                perturbed_accum_error = self.AccumulatedError_J(y, y_hat_perturbed)
                self.w1[i,j]         -= perturbation
                self.dJ_dw1[i,j]      = (perturbed_accum_error - current_accumulated_error) / perturbation

        # compute gradients of the wij^(2)
        rows, cols = neural_network.w2.shape
        for i in range(rows):
            for j in range(cols):
                self.w2[i,j]         += perturbation
                y_hat_perturbed       = self.Forward(x)
                perturbed_accum_error = self.AccumulatedError_J(y, y_hat_perturbed)
                self.w2[i,j]         -= perturbation
                self.dJ_dw2[i,j]      = (perturbed_accum_error - current_accumulated_error) / perturbation
    
    def UpdateWeightsSteepestDescent(self, steepest_descent_factor):
        self.w1 -= self.dJ_dw1 * steepest_descent_factor
        self.w2 -= self.dJ_dw2 * steepest_descent_factor


# We define the raw data: y=x**2
x = np.transpose(np.array(([np.linspace(-10,10,50)]), dtype=float)) # needs to be an array with 2 dimensions
y = x # x**2

# We scale the data...
x_max = x.max()
y_max = y.max()
x = x / x_max
y = y / y_max
# plot reference data
pl.plot(x, y, "b +", label = "Data")

# now we plot the first estimation...
neural_network = NeuralNetwork()
y_hat = neural_network.Forward(x)
pl.plot(x, y_hat, color="k", linewidth=5, label = "Initial guess")

current_accumulated_error = neural_network.AccumulatedError_J(y, y_hat)

'''
General parameters
'''
echo_level = 0

'''
Numerical alg parameters
'''
perturbation            = 1.0e-9  # for computing the grandients
steepest_descent_factor = 1.0e-4 # For updating the weights

'''
Convergence parameters
'''
tolerance = 1.0e-8
max_iter  = 5000e3

# we initialize some values
iteration = 0
old_accumulated_error = 1.0
relative_error = 1.0

# interval of print info stream and plot
interval_iter_print = 2000 # iterations
print_counter = 0

while relative_error > tolerance and iteration < max_iter:
    iteration += 1
    # compute gradients of the wij^(1)
    neural_network.Compute_dJdw_finite_differences(perturbation, current_accumulated_error)

    # now we update the weights ==> Learning...
    neural_network.UpdateWeightsSteepestDescent(steepest_descent_factor)

    # new estimation
    y_hat = neural_network.Forward(x)

    current_accumulated_error = neural_network.AccumulatedError_J(y, y_hat)
    relative_error = np.abs((old_accumulated_error - current_accumulated_error) / current_accumulated_error)
    old_accumulated_error = current_accumulated_error

    print_counter += 1
    if print_counter > interval_iter_print:
        print(" ## Iteration: ", iteration)
        print("    The current relative_error is:    ", "{0:.4e}".format(relative_error).rjust(11), "\n")
        # print("    The current_accumulated_error is: ", "{0:.4e}".format(current_accumulated_error).rjust(11), "\n")
        print_counter = 0
        pl.plot(x, y_hat, color="r")

if iteration == max_iter:
    print("The algorithm did not converge in the max iterations, rel error achieved: ", "{0:.4e}".format(relative_error).rjust(11))

if echo_level > 0:
    neural_network.PrintWeights()

pl.xlabel('x', fontsize = 12)
pl.ylabel('y', fontsize = 12)
pl.plot(x, y_hat, color="g", linewidth=5, label="Final result")
pl.legend(fontsize = 13)
pl.grid()
pl.show()
