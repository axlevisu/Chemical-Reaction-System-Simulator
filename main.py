# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from fractions import Fraction

def arraypow(x,A):
	return np.prod(x**(A.T),axis=-1)

def ode(y,t,A,Ok):
	# theta, X = np.array(y[:A.shape[0]]), np.array(y[A.shape[0]:])
	theta,X = y[:A.shape[0]],y[A.shape[0]:]
	# theta,X = y
	# print y
	MprojReaction = A.dot(X - arraypow(theta,A))
	forward_rate = arraypow(theta,-1*A.dot(Ok*(Ok <0)))
	backward_rate = arraypow(theta,-1*A.dot(Ok*(Ok >0)))
	EProjReaction = Ok.dot(backward_rate*(arraypow(X,-Ok*(Ok <0))) -  forward_rate*arraypow(X,Ok*(Ok >0))) 
	# print MprojReaction.tolist() + EProjReaction.tolist()
	# return MprojReaction.tolist() + EProjReaction.tolist()
	# print np.concatenate((MprojReaction,EProjReaction))
	return np.concatenate((MprojReaction,EProjReaction)).tolist()

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[0,2,3],[1,1,1]]
u = [1,1]
X_init = [0.5,0.25,0.25]
param_init = [1,1]
# A = [[3,1,0,2],[0,2,3,1]]
# O = [[1,1,0,0],[1,1,1,1]]
# u = [1,1]
# X_init = [0.3,0.2,0.1,0.4]
# param_init = [1,1]

# Number of timesteps
ts = 10000

A = np.array(A)
O = np.array(O)
u = np.array(u)
Ker = 1.0*np.array(Matrix(O).nullspace())
# Making the entries in the kernel basis integers
Ok=[]
for basis in Ker: 
	l = lcm(map(lambda x: Fraction(x).denominator,map(str,basis)))
	basis = map(int,l*basis)
	Ok.append(basis)

# Kernel basis are columns
Ok = np.array(Ok).T
ts = ts
param_init = np.array(param_init)
X_init = np.array(X_init)
# Normalizing the parameters
if param_init is None:
	theta = 1.0*np.ones(A.shape[0])/A.shape[1]
else:
	# Making sure sum of theta^A is 1 initially
	theta = 1.0*param_init/np.sum(arraypow(param_init,A))

if X_init is None:
	# This code is incomplete
	# TODO: Add code here for an appropiate kernel element, use sympy for null space
	# seed = np.random.normal(0,1)
	B = np.linalg.pinv(O)
	X =  B.dot(u)
	# k = np.random.uniform(0,-np.min(X0))
else:
	# Makes sure sum of X is 1
	X = 1.0*X_init/sum(X_init)

print "Design Matrix:"
print A

print "Observation Matrix:"
print O

print "Kernel Basis of O:"
print Ok

print "initial theta and X:"
print theta, X

t = np.linspace(0, 100, ts)
# y0 = theta.tolist() + X.tolist()
y0 = np.concatenate((theta,X))
# y0 = (theta,X)
sol = odeint(ode, y0, t, args=(A,Ok))
# Final theta and X
y = sol[-1,:]
theta = sol[-1,:A.shape[0]]
X = sol[-1,A.shape[0]:]
print "Final Theta and  X:"
print theta,X
# Print the rates to check
print "Final derivatives:"
print ode(y,t,A,Ok)

