import numpy as np 


#This class just makes/initialises the nodes and the matrix A, O for a boltzmann machine so that we can run it 
#on our CRN simulator
class RBM:
	def __init__(self,visible_layer=2,hidden_layer=1,fully_connected=True):
		self.visible_layer = visible_layer
		self.hidden_layer = hidden_layer 
		self.fully_connected = fully_connected #fraction of nodes not connected

	def create_boltzmann_machineA(self):
		nos_independent_vars=self.visible_layer+self.hidden_layer
		nos_dependent_vars   = self.visible_layer*self.hidden_layer # = theta_{ij}
		total_theta = nos_independent_vars + nos_dependent_vars
		if (self.fully_connected):
			A = np.zeros((nos_independent_vars+nos_dependent_vars,2**nos_independent_vars),dtype=np.int8)
			for i in range(2**nos_independent_vars):
				rbmstate = bin(i)[2:]
				rbmstate = '0'*(nos_independent_vars-len(rbmstate)) + rbmstate  
				for j in range(nos_independent_vars+nos_dependent_vars):
					if (j < nos_independent_vars):
						A[j,i] = int(rbmstate[-(j+1)]);
					else : 
						quo = (j-nos_independent_vars) / self.hidden_layer  #independent
						rem = (j-nos_independent_vars) - (self.hidden_layer*quo) 
						A[j,i] = int(rbmstate[-(quo+1)])*int(rbmstate[-(self.visible_layer+rem+1)])
						# Making Column sums equal
 			last_row = (nos_independent_vars+nos_dependent_vars)*np.ones(2**nos_independent_vars,dtype=np.int8)
 			last_row = last_row - np.sum(A,axis=0)
 			A = np.r_[A,[last_row]]
		return A

	def create_boltzmann_machineO(self):
		O = np.zeros((2**self.visible_layer,2**(self.visible_layer+self.hidden_layer)),dtype=np.int8)
		for i in range(2**self.visible_layer):
			rbmstate1 = bin(i)[2:]
			rbmstate1 = '0'*(self.visible_layer - len(rbmstate1)) + rbmstate1
			for j in range(2**self.hidden_layer):
				rbmstate2 = bin(j)[2:]
				rbmstate2 = '0'*(self.hidden_layer - len(rbmstate2)) + rbmstate2
				rbmstate = rbmstate2 + rbmstate1
				int_rbmstate = int(rbmstate,2)
				O[i,int_rbmstate] = 1
		return O


def main():
	trivial = RBM()
	A = trivial.create_boltzmann_machineA()
	O = trivial.create_boltzmann_machineO()
	print A
	print O


if __name__ == '__main__':
	main()








