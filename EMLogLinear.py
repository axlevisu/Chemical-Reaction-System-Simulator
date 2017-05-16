# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def reduce_int(alist):
	mul = 1
	# abs_list = abs(alist)
	for x in abs_list:
		x = float(x - floor(x))
		if x:
			mul *= Fraction(x).denominator
	abs_list = [a*mul for a in abs_list]
	gcd = reduce(lambda x,y: gcd([x,y]),abs_list)
	mul = mul/gcd
	alist = [1*mul for a in alist]
	return alist

def reaction_eq(conc, coeff):
	eq = 1
	for pair in transpose(conc):
		eq *= float(pair[0])**pair[1]
	return eq


def transpose(alist):
	return map(list, zip(*alist))

def arraypow(x,A):
	return np.prod(x**(A.T),axis=-1)

class EMLogLinear():
	"""docstring for EMLogLinear"""
	def __init__(self, A,O,u,param_init = None, X_init = None,ts = 1000):
		# super(EMLogLinear, self).__init__()
		self.A = np.array(A)
		self.O = np.array(O)
		self.u = np.array(u)
		self.Ok = np.array(Matrix(O).nullspace())
		self.ts = ts
		param_init = np.array(param_init)
		X_init = np.array(X_init)
		# Normalizing the parameters
		if param_init is None:
			self.theta = 1.0*np.ones(self.A.shape[0])/self.A.shape[1]
		else:
			self.theta = 1.0*param_init/np.sum(arraypow(param_init,self.A))
		if X_init is None:
			# This code is incomplete
			# TODO: Add code here for an appropiate kernel element, use sympy for null space
			# seed = np.random.normal(0,1)
			B = np.linalg.pinv(self.O)
			self.X =  B.dot(self.u)
			# k = np.random.uniform(0,-np.min(self.X0))
		else:
			self.X = 1.0*X_init/sum(X_init)

	def ode(y,t):
		theta, X = y 
		MprojReaction = self.A.dot(X - arraypow(theta,self.A))
		rate = arraypow(theta,self.A.dot(self.Ok))
		EProjReaction = -self.Ok.dot( arraypow(X,self.Ok*(self.Ok >0)) - rate*(arraypow(X,self.Ok*(self.Ok <0))) ) 
		return [MprojReaction,EProjReaction]

	def solve(self):
		t = np.linspace(0, 50, self.ts)
		y0 = [self.theta, self.theta]
		sol = odeint(ode, y0, t)
		return sol[-1,:]
	# def EProj():
	# 	X0 = self.X
	# 	q = self.theta**(self.A.T)
		

def main():
	A = [[2,1,0],[0,1,2]]
	O = [[0,1,3],[1,1,1]]
	u = [1,1]
	X_init = [0.5,0.25,0.25]
	param_init = [1,1]
	model = EMLogLinear(A = A, O =O, u = u, X_init = X_init, param_init = param_init)
	print model.solve()
	return


if __name__ == '__main__':
	main()