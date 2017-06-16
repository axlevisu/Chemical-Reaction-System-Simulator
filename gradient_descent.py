# gradient_descent.py
# partial_loglinear.py
import numpy as np
from sympy import Matrix,lcm
import matplotlib.pyplot as plt
from fractions import Fraction
from timeit import default_timer

random_seed =6
np.random.seed(random_seed)

def Lyapunov(X,Y):
	return np.sum(Y - X + (X*np.log(X/Y)))

def KerIntBasis(B):
	BKer = 1.0*np.array(Matrix(B).nullspace())
	Bk =[]
	for basis in BKer: 
		l = lcm(map(lambda x: Fraction(x).limit_denominator().denominator,map(str,basis)))
		basis = map(int,l*basis)
		Bk.append(basis)	
	Bk = np.array(Bk).T #Basis are column elements
	return Bk

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def partial_gradient(theta,X,A,Ok):
	Y = arraypow(theta,A)
	Z = np.sum(Y)
	Y = Y/Z
	gt = -(A.dot(X - Y))/theta
	gX = (Ok.dot(1. +np.log((X)/Y))).dot(Ok)
	return gt,gX

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[7,13,5],[1,1,1]]
X_init = [0.1,0.7,0.2]
# A = [[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1],[0,0,0,0,0,0,1,1],[0,0,0,1,0,0,0,1],[5,4,4,2,4,3,2,0]]
# O = [[1,0,1,0,0,0,0,0],[0,1,0,1,0,0,0,0],[0,0,0,0,1,0,1,0],[0,0,0,0,0,1,0,1]]
# X_init = [1/6.,1/6.,1/6.,1/6.,1/12.,1/12.,1/12.,1/12.]
# param_init = [1.,1.,1.,1/.2,1/.2,1.]
A = np.array(A)
O = np.array(O)
# Randomly initialize parameters
theta = np.random.uniform(0.01,1,A.shape[0])
X = np.array(X_init)/np.sum(X_init)
u = O.dot(X)
Ok = KerIntBasis(O).T
Ok = Ok*2.0/np.max(np.abs(Ok)) 
Ak = KerIntBasis(A).T


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
al= 0.0001
# Niter = 100000
lamda =1
eps = 10**(-9)
print "Learning Rate:",al
print "Gradient treshold:", eps
Niter =0
start = default_timer()
max_iter =100000
# while(True):
for i in xrange(max_iter):
	gt,gX = partial_gradient(theta,X,A,Ok)
	if (np.abs(gt)<eps).all() and (np.abs(gX)<eps).all():
		break
	theta = theta -al*gt
	X = X-al*gX
	Niter +=1
stop = default_timer()
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

print "Gradient Descent Took:", str(stop-start) +"s"