# partial_loglinear.py
import numpy as np
from sympy import Matrix,lcm
import matplotlib.pyplot as plt
from fractions import Fraction
from CRN.ReactionSystem import MassActionSystem, arraypow
from timeit import default_timer
from RBM.RBM import RBM
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

# Give Model Here
A = [[2,1,0],[0,1,2]]
O = [[7,10,2],[1,1,1]]
X_init = [0.1,0.7,0.2]
# A = [[3,1,0,2],[0,2,3,1]]
# O = [[1,1,0,0],[1,1,1,1]]
# X_init = [0.3,0.2,0.1,0.4]
# A = [[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1],[0,0,0,0,0,0,1,1],[0,0,0,1,0,0,0,1],[5,4,4,2,4,3,2,0]]
# O = [[1,0,1,0,0,0,0,0],[0,1,0,1,0,0,0,0],[0,0,0,0,1,0,1,0],[0,0,0,0,0,1,0,1]]
# X_init = [1/6.,1/6.,1/6.,1/6.,1/12.,1/12.,1/12.,1/12.]

# Number of timesteps
ts = 1000000
t = 100000

A = np.array(A)
O = np.array(O)
X_init = np.array(X_init)
Ok = KerIntBasis(O).T
Ok = Ok*2.0/np.max(np.abs(Ok)) 
Ak = KerIntBasis(A).T

# Randomly initialize parameters
theta = np.random.uniform(0.01,1,A.shape[0])
# theta = 1.0*np.ones(A.shape[0])

if X_init is None:
	B = np.linalg.pinv(O)
	X =  B.dot(u) # TODO: Add a vector from the nullspace and make sure X is positive
else:
	# Makes sure sum of X is 1
	X = 1.0*X_init/sum(X_init)

u = O.dot(X)
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
species =[]
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

for i in xrange(A.shape[0]):
	species.append("t_"+str(i+1))

for i in xrange(O.shape[1]):
	species.append("X_"+str(i+1))

reactions = np.array(reactions)
rates = np.array(rates)
system = MassActionSystem(reactions,rates,species)
print "The Reactions are:"
print system.display_reactions()
system.set_concentrations(Y_init)
# output = system.run(t=t,ts=ts)
delta_time =0.001
eps = 10**-12
print "Each time step is:",str(delta_time)+"s"
print "Gradient treshold:", eps
start = default_timer()
# output,t = system.run_till(delta_time=delta_time,eps=eps,every=100)
output = system.run(t=1000,ts=100000)
# output = system.run()
stop = default_timer()
Y = system.current_concentrations()
# print "Ran for:", str(t)+"s", "or",t/delta_time,"iterations"
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
print "Mass Action Kinetics Took:", str(stop-start) +"s"
# t = np.linspace(0, t, ts)
# lyapunov = [Lyapunov(arraypow(y[:St],A),y[St:]) for y in output]
# plt.plot(t, lyapunov, label="Lyapunov Function")
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.grid()
# plt.show()

