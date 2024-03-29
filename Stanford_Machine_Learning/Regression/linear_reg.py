
import numpy as np
import matplotlib.pyplot as pl

"""  Supervised learning!
In this file we create a linear regression algorithm to model a random data set
"""

# --------------------------------------------------------------
def h(theta1, theta2, x):
    return theta1 + x * theta2

# --------------------------------------------------------------
def Calculate_J_error(theta1, theta2, x_i, y_i):
    square_sum = 0.0
    for i in range(x_i.size):
        square_sum += (h(theta1, theta2, x_i[i]) - y_i[i])**2
    return 0.5 * square_sum

# --------------------------------------------------------------
def Calculate_dJ_dTheta(theta1, theta2, x_i, y_i):
    dJ_dTheta = np.zeros(2) # dJ_dT1 , dJ_dT2

    for i in range(x_i.size):
        factor = h(theta1, theta2, x_i[i]) - y_i[i]
        dJ_dTheta[0] += factor
        dJ_dTheta[1] += factor * x_i[i]

    return dJ_dTheta
# --------------------------------------------------------------
# Create the data set

x = np.random.rand(100)
y = 3.0 * x + np.random.randn(100)

pl.plot(x, y, "bo", label="dataset")


# initial guess
Theta1 = 1.0
Theta2 = 1.0

J     = Calculate_J_error(Theta1, Theta2, x, y)
dJ_dT = Calculate_dJ_dTheta(Theta1, Theta2, x, y)

print("*******************************")
print("The initial error J is: ", J)
print("The initial dJ_dT is  : ", dJ_dT)

tolerance = 1.0e-8
learning_rate = 0.01
iter = 0
max_iter = 2000
Jold = 0.0

while abs(J-Jold) > tolerance:
    Theta1   -= learning_rate * dJ_dT[0]
    Theta2   -= learning_rate * dJ_dT[1]
    Jold      = J
    J         = Calculate_J_error(Theta1, Theta2, x, y)
    dJ_dT     = Calculate_dJ_dTheta(Theta1, Theta2, x, y)
    print("Iteration: ", iter, "Error |J-Jold|: ", '{:.2e}'.format(abs(J-Jold)))
    iter += 1

    if iter >= max_iter:
        break
print("*******************************")
print("The linear equation is -->  y = ", '{:.2e}'.format(Theta1), " + ", '{:.2e}'.format(Theta2), " * x")

pl.plot(x, h(Theta1, Theta2, x), color="red", label="Linear regression")

#-----------------------------------------
# try now with stochastic gradient descent
#-----------------------------------------
StoTheta1 = 0.0
StoTheta2 = 0.0
iter = 0
max_iter = 2500
Jold = 0.0
J    = Calculate_J_error(StoTheta1, StoTheta2, x, y)
learning_rate = 0.01
tolerance = 1.0e-6

print("*******************************")
print("The initial Stochastic error J is: ", J)

while abs(J-Jold) > tolerance:
    for i in range(x.size): # we update the parameters only according to the derivative of ONE sample at each time (no sumation)
        factor     = h(StoTheta1, StoTheta2, x[i]) - y[i]
        StoTheta1 -= learning_rate * factor * 1.0
        StoTheta2 -= learning_rate * factor * x[i]
        Jold      = J
        J         = Calculate_J_error(StoTheta1, StoTheta2, x, y)

        print("Stochastic Iteration: ", iter, "Error |J-Jold|: ", '{:.2e}'.format(abs(J-Jold)))
        iter += 1
        if iter >= max_iter:
            break
    if iter >= max_iter:
        break
print("*******************************")
print("The linear equation is -->  y = ", '{:.2e}'.format(StoTheta1), " + ", '{:.2e}'.format(StoTheta2), " * x")

pl.plot(x, h(StoTheta1, StoTheta2, x), color="g", label="Stochastic Linear regression")

pl.legend()
pl.show()