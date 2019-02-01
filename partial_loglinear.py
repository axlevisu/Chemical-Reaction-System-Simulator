# partial_loglinear.py
import numpy as np
from sympy import Matrix,lcm
import matplotlib.pyplot as plt
from fractions import Fraction
from CRN.ReactionSystem import MassActionSystem, arraypow
from timeit import default_timer
from RBM.RBM import RBM
from LLL import basis_reduce
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

# Number of timesteps
ts = 2000
t =  0.8

#Load from file
A =[]
O =[]
model_file_name = "model-observation.txt"
model_file = open(model_file_name, 'r')
St = int(model_file.readline().strip())
Sx = int(model_file.readline().strip())
for i in xrange(St):
	line = model_file.readline()
	A.append(list(map(int, line.strip().split(','))))

Sr = int(model_file.readline().strip())
for i in xrange(Sr):
	line = model_file.readline()
	O.append(list(map(int, line.strip().split(','))))

line = model_file.readline()
print line
X_init = list(map(float, line.strip().split(',')))

A = np.array(A)
O = np.array(O)
# X_init = 1.0*np.array(X_init)/np.sum(X_init)
X_init = 1.0*np.array(X_init)
# X = X_init/np.sum(X_init)
X = X_init
print "Before LLL"
Ok = KerIntBasis(O).T
Ak = KerIntBasis(A).T
print Ok
Ok =basis_reduce(Ok)
# Randomly initialize parameters
theta = np.random.uniform(0.01,1,A.shape[0])
theta = np.ones(A.shape[0])
# Calcumating u
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
# Ok = Ok*1./np.max(np.abs(Ok)) 
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
	species.append("\\theta_"+str(i+1))

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
output = system.run(t=t,ts=ts, plot=True)
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
