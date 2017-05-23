# partial_loglinear.py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from fractions import Fraction
from ReactionSystem.MassActionSystem import MassActionSystem, arraypow

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
ts = 1000
t = 10

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

reactions =[]
rates =[]
Y_init = np.concatenate((theta,X))
S = Y_init.shape[0] # Number of species
St = theta.shape[0] # Number of theta species
Sx = X.shape[0] # Number of X species
# 2 for every column
for i in xrange(A.shape[1]):
	column = A.T[i].tolist()
	complexl = column + [0]*Sx
	complexr = [0]*S
	reactions.append([complexl[:],complexr[:]])
	complexr[i+St] =1
	complexl[i+St] =1
	reactions.append([complexr,complexl])
	rates = rates +[1.,1.]

for i in xrange(Ok.shape[0]):
	kernel = Ok[i]
	catalysis = kernel.dot(A.T)
	f = (-catalysis*(catalysis <0)).tolist() #Forward reaction catalyst
	b = (catalysis*(catalysis >0)).tolist() #Backward reaction catalyst
	p = (kernel*(kernel>0)).tolist()                #Positive complex
	n = (-kernel*(kernel<0)).tolist()				#Negative complex
	reactions = reactions + [ [f + p,f + n],[b + n,b + p]]
	rates = rates +[1.,1.]

reactions = np.array(reactions)
rates = np.array(rates)
system = MassActionSystem(reactions,rates)
system.set_concentrations(Y_init)
output = system.run(t=t,ts=ts)
Y = system.current_concentrations()

theta = Y[:St]
X =Y[St:]
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

print "Final dy/dt:"
print system.dydt()

print "Final Lyapunov"
print Lyapunov(X,MLD)

t = np.linspace(0, t, ts)
lyapunov = [Lyapunov(arraypow(y[:St],A),y[St:]) for y in output]
plt.plot(t, lyapunov, label="Lyapunov Function")
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()