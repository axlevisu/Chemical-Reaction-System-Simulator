# compare.py
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from fractions import Fraction
from CRN.ReactionSystem import MassActionSystem, arraypow
from timeit import default_timer
from RBM.RBM import RBM
from LLL import basis_reduce
# random_seed =6
# np.random.seed(random_seed)

def multichoose(n,k):
    if k < 0 or n < 0: return "Error"
    if not k: return [[0]*n]
    if not n: return []
    if n == 1: return [[k]]
    return [[0]+val for val in multichoose(n-1,k)] + [[val[0]+1]+val[1:] for val in multichoose(n,k-1)]

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

def partial_gradient(theta,X,A,Ok):
	Y = arraypow(theta,A)
	Z = np.sum(Y)
	Y = Y/Z
	gt = -(A.dot(X - Y))/theta
	gX = (Ok.dot(1 + np.log((X)/Y))).dot(Ok)
	return gt,gX

# Load model
A =[]
O =[]
model_file_name = "model-observation-ABC.txt"
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


A = np.array(A)
O = np.array(O)
Ok = KerIntBasis(O).T
Ok = basis_reduce(Ok)
Ok = Ok*1./np.max(np.abs(Ok))  
Ak = KerIntBasis(A).T
reactions =[]
rates =[]
S = Sx+St
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
delta_time =0.001
eps = 10**-9
al= 0.0001
max_iter =400000
print "O:"
print O
print "A:"
print A
print "al:",al
print "max_iter:",max_iter
print "gradient_eps",eps
t =60000
ts=1000000
print "t:",t,"ts:",ts
# X_init = np.array(multichoose(Sx,20))/20.0
X_init = np.random.rand(50,Sx)
X_init = (X_init.T/np.sum(X_init,axis=-1)).T
# output_file_name = 
print X_init.shape
total=0
crn_better =0
# X_init = np.array([[0.9,0.025,0.075]])
# for X in np.flip(X_init,0):
for X in X_init:
	if not (X==0).any():
		theta = np.random.uniform(0.01,1,A.shape[0])
		# theta = np.array([ 0.17402282,  0.39597236])
		Y_init = np.concatenate((theta,X))
		system.set_concentrations(Y_init)
		# output,t = system.run_till(delta_time=delta_time,eps=eps,every=1000)
		output = system.run(t=t,ts=ts)
		Y = system.current_concentrations()
		crn_theta = Y[:St]
		crn_X =Y[St:]
		crn_MLD =arraypow(crn_theta,A)
		CRN_Lp = Lyapunov(crn_X,crn_MLD)
		Xi = X
		thetai = theta
		for i in xrange(max_iter):
			gt,gX = partial_gradient(theta,X,A,Ok)
			if (np.abs(gt)<eps).all() and (np.abs(gX)<eps).all():
				break
			theta = theta -al*gt
			X = X-al*gX
		# print gt,gX
		MLD =arraypow(theta,A)
		MLD = MLD/np.sum(MLD)
		gd_Lp = Lyapunov(X,MLD)
		diff = CRN_Lp-gd_Lp
		if diff==nan:
			diff = 0
		
		if diff <=0 or abs(diff)<eps:
			diff = True
		else:
			diff =False
		print thetai,",",Xi,",",CRN_Lp,",",gd_Lp,",",diff
		total+=1
		crn_better+=diff

print crn_better,total
print 100.0*crn_better/total
		# print X,crn_X
		# print MLD, crn_MLD
