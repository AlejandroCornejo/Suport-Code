
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
        self.hidden_layer_size = 5  # number of hidden neurons

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
    def SigmoidDerivative(self, z):
        return np.exp(-z) / ((1.0 + np.exp(-z))**2)

    def Forward(self, x):
        self.z2 = np.dot(x, self.w1)        # Activity:   number_of_examples x number_of_hidden_layers
        self.a2 = self.Sigmoid(self.z2)     # Activation: number_of_examples x number_of_hidden_layers
        self.z3 = np.dot(self.a2, self.w2)  # Activity: 
        self.y_hat = self.Sigmoid(self.z3)  # prediction
        return self.y_hat

    def AccumulatedError_J(self, y, y_hat):
        # let's print the accumulated error ==>  AccumErr = Sum(0.5*(y-yhat)^2)
        error_y = y - y_hat
        return 0.5*(np.linalg.norm(error_y)**2)
    
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


    def Compute_dJdw_analytical(self, x, y):
        # compute dJ/dw2
        y_hat = self.Forward(x)
        delta3 = np.multiply(-(y-y_hat), self.SigmoidDerivative((self.z3)))
        self.dJ_dw2 = np.dot(np.transpose(self.a2), delta3)

        # compute dJ/dw1
        delta2 = np.dot(delta3, np.transpose(self.w2))*self.SigmoidDerivative((self.z2))
        self.dJ_dw1 = np.dot(np.transpose(x), delta2)
    
    def UpdateWeightsGradientDescent(self, gradient_descent_factor):
        self.w1 -= self.dJ_dw1 * gradient_descent_factor
        self.w2 -= self.dJ_dw2 * gradient_descent_factor


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
# perturbation            = 1.0e-9  # for computing the grandients
gradient_descent_factor = 1.0 # For updating the weights

'''
Convergence parameters
'''
tolerance = 1.0e-2
max_iter  = 1e9

# we initialize some values
iteration = 0

# interval of print info stream and plot
interval_iter_print = 5000 # iterations
print_counter = 0
is_converged = False

while not is_converged and iteration < max_iter:
    iteration += 1
    # compute gradients of the wij^(1)
    # neural_network.Compute_dJdw_finite_differences(perturbation, current_accumulated_error)
    neural_network.Compute_dJdw_analytical(x, y)

    # now we update the weights ==> Learning...
    neural_network.UpdateWeightsGradientDescent(gradient_descent_factor)

    # new estimation
    y_hat = neural_network.Forward(x)

    current_accumulated_error = neural_network.AccumulatedError_J(y, y_hat)

    norm_derivatives = np.linalg.norm(neural_network.dJ_dw1) + np.linalg.norm(neural_network.dJ_dw2)
    is_converged = (norm_derivatives) <= tolerance

    print_counter += 1
    if print_counter > interval_iter_print:
        print(" ## Iteration: ", iteration)
        print("    The current norm_derivatives is:          ", "{0:.4e}".format(norm_derivatives).rjust(11))
        print("    The current current_accumulated_error is: ", "{0:.4e}".format(current_accumulated_error).rjust(11), "\n")
        print_counter = 0
        pl.plot(x, y_hat, color="r")

if iteration == max_iter:
    print("The algorithm did not converge in the max iterations, norm_derivatives achieved: ", "{0:.4e}".format(norm_derivatives).rjust(11))
else:
    print("The algorithm converged in ", str(iteration) ," iterations, norm_derivatives achieved: ", "{0:.4e}".format(norm_derivatives).rjust(11))

if echo_level > 0:
    neural_network.PrintWeights()
    neural_network.PrintGradients()

pl.xlabel('x', fontsize = 15)
pl.ylabel('y', fontsize = 15)
pl.plot(x, y_hat, color="g", linewidth=5, label="Final result")
pl.legend(fontsize = 13)
pl.grid()
pl.show()
