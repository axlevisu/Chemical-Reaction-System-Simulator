# EM-Log Linear
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def reduce_int(alist):
	mul = 1
	abs_list = abs(alist)
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


class EMLogLinear():
	"""docstring for EMLogLinear"""
	def __init__(self, A,O,u,param_init = None, X_init = None,ts = 0.001):
		# super(EMLogLinear, self).__init__()
		self.A = np.array(A)
		self.O = np.array(O)
		self.u = np.array(u)
		self.Ok = np.array(Matrix(O).nullspace())
		if param_init is None:
			self.theta = np.ones(self.A.shape[0])
		else:
			self.theta = param_init
		if X_init is None:
			# TODO: Add code here for a appropiate kernel element, use sympy for null space
			# seed = np.random.normal(0,1)
			B = np.linalg.pinv(self.O)
			self.X =  B.dot(self.u)
			# k = np.random.uniform(0,-np.min(self.X0))
		else:
			self.X = X_init				

	# def EProj():
	# 	X0 = self.X
	# 	q = self.theta**(self.A.T)
		

def main():
	A = [[2,1,0],[0,1,2]]
	O = [[2,0,3],[1,1,1]]
	u = [1.4,1]
	model = EMLogLinear(A = A,O =O,u = u)
	return


if __name__ == '__main__':
	main()