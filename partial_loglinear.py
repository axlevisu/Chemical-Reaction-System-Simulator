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

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[0,1,3],[1,1,1]]
u = [1,1]
X_init = [0.5,0.25,0.25]
param_init = [1.,1.]
# A = [[3,1,0,2],[0,2,3,1]]
# O = [[1,1,0,0],[1,1,1,1]]
# X_init = [0.3,0.2,0.1,0.4]
# A = [[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1],[0,0,0,0,0,0,1,1],[0,0,0,1,0,0,0,1],[5,4,4,2,4,3,2,0]]
# O = [[1,0,1,0,0,0,0,0],[0,1,0,1,0,0,0,0],[0,0,0,0,1,0,1,0],[0,0,0,0,0,1,0,1]]
# X_init = [1/.6,1/.6,1/.6,1/.6,1/.12,1/.12,1/.12,1/.12]
# param_init = [1.,1.,1.,1/.2,1/.2,1.]

# Number of timesteps
ts = 1000000
t = 100000

A = np.array(A)
O = np.array(O)
u = np.array(u)
param_init = np.array(param_init)
X_init = np.array(X_init)
Ok = KerIntBasis(O).T 
Ak = KerIntBasis(A).T

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
print O.dot(X)

reactions =[]
rates =[]
Y_init = np.concatenate(theta,X)
S = Y_init.shape[0] # Number of species
St = theta.shape[0] # Number of theta species
Sx = X.shape[0] # Number of X species
# 2 for every column
for i in xrange(A.shape[1]):
	column = A.T[i]
	complexl, complexr = np.zeros(S),np.zeros(S)
	complexl[0:St] = column
	reactions.append([complexl,complexr])
	complexr[i+Sx] =1
	complexl[i+Sx] =1
	reactions.append([complexr,complexl])
	rates = rates +[1.,1.]

for i in xrange(Ok.shape[0]):
	kernel = Ok[i]
	catalysis = kernel.dot(A.T)
	complexl, complexr = np.zeros(S),np.zeros(S)
	complexl[St:] =  Ok*(Ok>0) 
	complexr[St:] = -Ok*(Ok<0)
	complexl[:St] = -catalysis*(catalysis <0)
	complexr[:St] = -catalysis*(catalysis <0)
	reactions.append([complexl,complexr])
	complexl[:St] = catalysis*(catalysis >0)
	complexr[:St] = catalysis*(catalysis >0)
	reactions.append([complexr,complexl])
	rates = rates +[1.,1.]

reactions = np.array(reactions)
rates = np.array(rates)
Y_init = np.concatenate(theta,X)

