# EMmain.py
# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from fractions import Fraction

def arraypow(x,A):
	return np.prod(x**(A.T),axis=-1)

def EProj(y,t,Ok):
	Y,X = y[:Ok.shape[0]],y[Ok.shape[0]:]
	forward_rate = arraypow(Y,-Ok*(Ok <0))
	backward_rate = arraypow(Y,Ok*(Ok >0))
	EProjReaction = Ok.dot(backward_rate*(arraypow(X,-Ok*(Ok <0))) -  forward_rate*arraypow(X,Ok*(Ok >0))) 
	# print EProjReaction
	return [0.]*Ok.shape[0] + EProjReaction.tolist() 

def MProj(y,t,Ak):
	Y,X = y[:Ak.shape[0]],y[Ak.shape[0]:]
	MprojReaction = Ak.dot( arraypow(X,-Ak*(Ak <0)) -  arraypow(X,Ak*(Ak >0)) )
	# print MprojReaction
	return [0.]*Ak.shape[0] + MprojReaction.tolist() 

def KLDiv(Y,X):
	return np.sum(X*np.log(X/Y))
# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[0,2,3],[1,1,1]]
u = [1.25,1]
X_init = [0.5,0.25,0.25]
param_init = [1,1]
# A = [[3,1,0,2],[0,2,3,1]]
# O = [[1,1,0,0],[1,1,1,1]]
# X_init = [0.3,0.2,0.1,0.4]

# Number of timesteps
ts = 10000

A = np.array(A)
O = np.array(O)
u = np.array(u)
OKer = 1.0*np.array(Matrix(O).nullspace())
AKer = 1.0*np.array(Matrix(A).nullspace())
# Making the entries in the kernel basis integers
Ok=[]
Ak=[]
for basis in OKer: 
	l = lcm(map(lambda x: Fraction(x).denominator,map(str,basis)))
	basis = map(int,l*basis)
	Ok.append(basis)

for basis in AKer: 
	l = lcm(map(lambda x: Fraction(x).denominator,map(str,basis)))
	basis = map(int,l*basis)
	Ak.append(basis)

# Kernel basis are columns
Ok = np.array(Ok).T
Ak = np.array(Ak).T
ts = ts
param_init = np.array(param_init)
X_init = np.array(X_init)
# Normalizing the parameters
if param_init is None:
	theta = 1.0*np.ones(A.shape[0])/A.shape[1]
else:
	# NOT NEEDED: Making sure sum of theta^A is 1 initially
	theta = 1.0*param_init

if X_init is None:
	# This code is incomplete
	# TODO: Add code here for an appropiate kernel element, use sympy for null space
	B = np.linalg.pinv(O)
	X =  B.dot(u) # TODO: Add a vector from the nullspace and make sure X is positive
else:
	# Makes sure sum of X is 1
	X = 1.0*X_init/sum(X_init)

print "Design Matrix:"
print A

print "Observation Matrix:"
print O

print "Kernel Basis of O:"
print Ok

print "Kernel Basis of A:"
print Ak


t = np.linspace(0, 100, ts)
Y = arraypow(theta,A)
Y = Y/np.sum(Y)
print "Starting point on toric ray:"
print Y
print "Starting point on affine space:"
print X
y0 = np.concatenate((Y,X))
eps = 10**-15
while(True):
	sol = odeint(EProj, y0, t, args=(Ok,))
	y = sol[-1,:]
	# if (y == y0).all():
	if KLDiv(y,y0)<eps:
		break
	# print "After EProjection:"
	# print y,y0
	y0 = y	
	sol = odeint(MProj, y0, t,args = (Ak,))
	y = np.concatenate((sol[-1,A.shape[1]:],y0[O.shape[1]:]))

	# print "After MProjection:"
	# print y,y0
	# if (y == y0).all():
	if KLDiv(y,y0)<eps:
		break
	y0 = y

# Final theta and X
y = sol[-1,:]
Y = sol[-1,:A.shape[1]]
X = sol[-1,O.shape[1]:]
print "Final MLD and  X:"
print Y,X

# Calculating KL - Divergence
print "KL Divergence:"
print KLDiv(Y,X)

print (Y[1]**2)/(Y[0]*Y[2])
print (X[1]**2)/(X[0]*X[2])
print O.dot(X)
print O.dot(Y)
