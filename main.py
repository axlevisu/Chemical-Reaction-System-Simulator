#Code for solving ode for  Chemical Reaction networks
# Here we are trying to solve for the following CRN
# We have the following reaction network

# The second order differential equation for the angle `theta` of a
# pendulum acted on by gravity with friction can be written::

# theta''(t) + b*theta'(t) + c*sin(theta(t)) = 0

# where `b` and `c` are positive constants, and a prime (') denotes a
# derivative.  To solve this equation with `odeint`, we must first convert
# it to a system of first order equations.  By defining the angular
# velocity ``omega(t) = theta'(t)``, we obtain the system::

# theta'(t) = omega(t)
# omega'(t) = -b*omega(t) - c*sin(theta(t))

# Let `y` be the vector [`theta`, `omega`].  We implement this system
# in python as:

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
def pend(y, t):
    x1,x2,x3,theta1, theta2 = y
    t1 = 2*(-theta1*(x1**2)*x3 + theta2*x2**3)
    t2 = -3*(-theta1*(x1**2)*x3 + theta2*x2**3)
    t3 = 1*(-theta1*(x1**2)*x3 + theta2*x2**3)
    t4 = 2*(x1-(theta1**2)) + (x2 - theta1*theta2)
    t5 = 2*(x3-(theta2**2)) + (x2 - theta1*theta2)
    dydt = [t1,t2,t3,t4,t5]
    return dydt
# ...

# We assume the constants are `b` = 0.25 and `c` = 5.0:
# For initial conditions, we assume the pendulum is nearly vertical
# with `theta(0)` = `pi` - 0.1, and it initially at rest, so
# `omega(0)` = 0.  Then the vector of initial conditions is

y0 = [0.5, 0.25,0.25,1,1]

# We generate a solution 101 evenly spaced samples in the interval
# 0 <= `t` <= 10.  So our array of times is:

t = np.linspace(0, 1, 1000)

# Call `odeint` to generate the solution.  To pass the parameters
# `b` and `c` to `pend`, we give them to `odeint` using the `args`
# argument.

sol = odeint(pend, y0, t)
plt.plot(t, sol[:, 0], 'b', label='x1')
plt.plot(t, sol[:, 1], 'g', label='x2')
plt.plot(t, sol[:, 2], 'r', label='x3')
plt.plot(t, sol[:, 3], 'c', label='\theta1')
plt.plot(t, sol[:, 4], 'm', label='\theta2')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()