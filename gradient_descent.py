# gradient_descent.py
# partial_loglinear.py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from fractions import Fraction

def Lyapunov(X,Y):
	return np.sum(Y - X + (X*np.log(X/Y)))

def KerIntBasis(B):
	BKer = 1.0*np.array(Matrix(B).nullspace())
	Bk =[]
	for basis in BKer: 
		l = lcm(map(lambda x: Fraction(x).denominator,map(str,basis)))
		basis = map(int,l*basis)
		Bk.append(basis)	
	Bk = np.array(Bk).T #Basis are column elements
	return Bk

def partial_gradient()

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[0,1,3],[1,1,1]]
X_init = [0.5,0.25,0.25]
param_init = [1.,1.]
A = np.array(A)
O = np.array(O)
param_init = np.array(param_init)
X_init = np.array(X_init)
u = O.dot(X_init)
Ok = KerIntBasis(O) 
Ak = KerIntBasis(A)

# Normalizing the parameters
if param_init is None:
	theta = 1.0*np.ones(A.shape[0])/A.shape[1]
else:
	# NOT NEEDED: Making sure sum of theta^A is 1 initially
	theta = 1.0*param_init

if X_init is None:
	B = np.linalg.pinv(O)
	X =  B.dot(u) # TODO: Add a vector from the nullspace and make sure X is positive
else:
	# Makes sure sum of X is 1
	X = 1.0*X_init/sum(X_init)

print "Design Matrix:"
print A

print "Observation Matrix:"
print O

print "Kernel Basis of O:(transposed)"
print Ok

print "Kernel Basis of A:(transposed)"
print Ak

print "initial theta and X:"
print theta, X

print "u (equals OX_init):"
print u

# Gradient descent params
al= 0.001
