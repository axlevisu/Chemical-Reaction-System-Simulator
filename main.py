# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from fractions import Fraction

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def ode(y,t,A,Ok):
	theta,X = y[:A.shape[0]],y[A.shape[0]:]
	MprojReaction = A.dot(X - arraypow(theta,A))
	forward_rate = arraypow(theta,-1*A.dot(Ok*(Ok <0)))
	backward_rate = arraypow(theta,1*A.dot(Ok*(Ok >0)))
	EProjReaction = Ok.dot(backward_rate*(arraypow(X,-Ok*(Ok <0))) -  forward_rate*arraypow(X,Ok*(Ok >0))) 
	return np.concatenate((MprojReaction,EProjReaction)).tolist()

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[0,1,3],[1,1,1]]
u = [1,1]
X_init = [0.6,0.1,0.3]
param_init = [0.5,0.2]
# A = [[3,1,0,2],[0,2,3,1]]
# O = [[1,1,0,0],[1,1,1,1]]
# X_init = [0.3,0.2,0.1,0.4]
# A = [[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1],[0,0,0,0,0,0,1,1],[0,0,0,1,0,0,0,1],[5,4,4,2,4,3,2,0]]
# O = [[1,0,1,0,0,0,0,0],[0,1,0,1,0,0,0,0],[0,0,0,0,1,0,1,0],[0,0,0,0,0,1,0,1]]
# u = [1/.3,1/.3,1/.6,1/.6]
# X_init = [1/.6,1/.6,1/.6,1/.6,1/.12,1/.12,1/.12,1/.12]
# param_init = [1.,1.,1.,1/.2,1/.2,1.]
# Number of timesteps
ts = 1000000

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
param_init = np.array(param_init)
X_init = np.array(X_init)
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

print "Kernel Basis of O:"
print Ok

print "Kernel Basis of A:"
print Ak

print "initial theta and X:"
print theta, X

print "u (equals OX_init):"
print O.dot(X)

t = np.linspace(0, 100000, ts)
y0 = np.concatenate((theta,X))
sol = odeint(ode, y0, t, args=(A,Ok))
# Final theta and X
y = sol[-1,:]
theta = sol[-1,:A.shape[0]]
X = sol[-1,A.shape[0]:]
print "Final Theta and  X:"
print theta,X

print "Y (MLD):"
Y = arraypow(theta,A)
print Y
# Print the rates to check
print "Final derivatives:"
print ode(y,t,A,Ok)
# Calculating KL - Divergence
print "KL Divergence:"
print np.sum(X*np.log(X/Y))
print "OX (should be equal to u):"
print O.dot(X)

print "OY, where Y is MLD"
print O.dot(Y)

print "X^Ak"
print arraypow(X,Ak)

params = A.shape[0];
for i in range(sol.shape[1]):
	if (i < params):
		labeli = 'theta' + str(i);
		plt.plot(t, sol[:, i], label=labeli)
	else :
		labeli = 'X' + str(i+1-params)
		plt.plot(t, sol[:, i],label=labeli)
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
