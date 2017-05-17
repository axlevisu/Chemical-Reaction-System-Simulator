# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from fractions import Fraction

# TODO: Possible to make tham member functions?
def ode(y,t,A,Ok):
	theta,X = y[:A.shape[0]],y[A.shape[0]:]
	MprojReaction = A.dot(X - arraypow(theta,A))
	forward_rate = arraypow(theta,-1*A.dot(Ok*(Ok <0)))
	backward_rate = arraypow(theta,-1*A.dot(Ok*(Ok >0)))
	EProjReaction = Ok.dot(backward_rate*(arraypow(X,-Ok*(Ok <0))) -  forward_rate*arraypow(X,Ok*(Ok >0))) 
	return np.concatenate((MprojReaction,EProjReaction)).tolist()

def arraypow(x,A):
	return np.prod(x**(A.T),axis=-1)

class EMLogLinear():
	"""docstring for EMLogLinear"""
	def __init__(self, A,O,u,param_init = None, X_init = None,ts = 1000):
		# super(EMLogLinear, self).__init__()
		self.A = np.array(A)
		self.O = np.array(O)
		self.u = np.array(u)
		Ker = 1.0*np.array(Matrix(O).nullspace())
		# Making the entries in the kernel basis integers
		self.Ok=[]
		for basis in Ker: 
			l = lcm(map(lambda x: Fraction(x).denominator,map(str,basis)))
			basis = map(int,l*basis)
			self.Ok.append(basis)

		# Kernel basis are columns
		self.Ok = np.array(self.Ok).T
		self.ts = ts
		param_init = np.array(param_init)
		X_init = np.array(X_init)
		# Normalizing the parameters
		if param_init is None:
			self.theta = 1.0*np.ones(self.A.shape[0])/self.A.shape[1]
		else:
			# NOT NEEDED: Making sure sum of theta^A is 1 initially
			self.theta = 1.0*param_init		
		if X_init is None:
			# This code is incomplete
			# TODO: Add code here for an appropiate kernel element, use sympy for null space
			B = np.linalg.pinv(self.O)
			self.X =  B.dot(self.u) # TODO: Add a vector from the nullspace and make sure X is positive
		else:
			# Makes sure sum of X is 1
			self.X = 1.0*X_init/sum(X_init)

	def solve(self):
		t = np.linspace(0, 50, self.ts)
		y0 = np.concatenate((self.theta,self.X))
		sol = odeint(ode, y0, t, args=(self.A,self.Ok))
		self.theta = np.array(sol[-1,:self.A.shape[0]])
		self.X = sol[-1,self.A.shape[0]:]
		return self.theta, self.X

	def divergence(self):
		Y = arraypow(self.theta,self.A)
		return np.sum(X*np.log(X/Y))

	def observation_kernel(self):
		return self.Ok

	
def main():
	A = [[2,1,0],[0,1,2]]
	O = [[0,2,3],[1,1,1]]
	u = [1,1]
	X_init = [0.5,0.25,0.25]
	param_init = [1,1]
	model = EMLogLinear(A = A, O =O, u = u, X_init = X_init, param_init = param_init)
	final_theta, final_X = model.solve()
	print final_theta
	print final_X
	return


if __name__ == '__main__':
	main()