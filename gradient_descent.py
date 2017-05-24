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

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def partial_gradient(theta,X,A,Ok):
	Y = arraypow(theta,A)
	gt = -(A.dot(X - Y))/theta
	gX = (Ok.dot(np.log(X/Y))).dot(Ok)
	return gt,gX

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
print u

# Gradient descent params
al= 0.0005
# Niter = 100000
lamda =1
eps = 10**(-12)
print "Learning Rate:",al
print "Gradient treshold:", eps
Niter =0
while(True):
	gt,gX = partial_gradient(theta,X,A,Ok)
	if (np.abs(gt)<eps).all() and (np.abs(gX)<eps).all():
		break
	theta = theta -al*gt
	X = X-al*gX
	Niter +=1

print gt,gX
# print theta,X
print "Number of iteration taken:",Niter
MLD =arraypow(theta,A)

print "Final Theta and X:"
print theta, X

print "MLD: theta^A"
print  MLD

print "Final sum of MLD"
print np.sum(MLD)

print "Final OX, should be equal to", u
print O.dot(X)

print "Final X^Ak:"
print arraypow(X,Ak.T)

print "Final Oxtheta^A or OxMLD, is it ",u,"?"
print O.dot(MLD) 

print "Final gradients:"
print gt,gX

print "Final Lyapunov"
print Lyapunov(X,MLD)
