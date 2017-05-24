# MassActionSystem.py
import numpy as np
from scipy.integrate import odeint

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def odes(y,t,R,k):
	l = R[:,0]
	r = R[:,1]
	return ((r-l).T).dot(arraypow(y,l.T)*k) 


class MassActionSystem(object):
	"""docstring for MassActionSystem"""
	def __init__(self, reactions, rates):
		# super(MassActionSystem, self).__init__()
		"""Each reaction in reactions is 
		a pair of complex (l,r) in l --> r"""
		reactions = np.array(reactions)
		rates = np.array(rates)*1.0
		self.reactions = reactions
		self.rates = rates

#### TODO: Add a method to add more reactions
#### TODO: Add a method to visualize

	def set_concentrations(self,concentrations):
		concentrations = 1.0*np.array(concentrations)
		if concentrations.shape[0] == self.reactions.shape[2]:
			self.concentrations = concentrations
		else:
			print "Wrong Initialization: Shapes of concentrations and reactions dont match ", concentrations.shape[0], "!=", self.reactions.shape[1]
		return self.concentrations

	def dydt(self):
		"""
		Returns rates of concentrations change at the current concentration
		"""
		return odes(self.concentrations,0,self.reactions,self.rates)

	def current_concentrations(self):
		"""
		Returns the concentrations after the latest run
		"""
		return self.concentrations

	def run(self,t=100000,ts=1000000):
		"""
		Run it for time t with ts number of time steps
		Outputs the concentration profile for the run
		default values are 100000 and 1000000
		"""	
		t = np.linspace(0, t, ts)
		y = self.concentrations
		output = odeint(odes, y, t, args= (self.reactions,self.rates))
		self.concentrations = output[-1,:]
		return output

	def run_till(self,time_step=0.001,eps = 10**-12):
		"""
		Runs till the gradients are less than eps
		returns time taken with array of concentrations
		default values of time_step and eps are 10**-3 and 10**-12 
		"""
		t = 0
		ts = time_step
		y = self.concentrations
		output = [y]
		while True:
			gy = odes(y,t,self.reactions,self.rates)
			if (np.abs(gy) <eps).all():
				break
			dy = gy*ts
			y += dy
			t+=ts
			output.append(y)

		output = np.array(output)
		self.concentrations =y
		return output,t
		
	

def main():
	reactions = [[[1,0],[0,1]], [[0,1],[1,0]]]
	rates = [1,1]
	system = MassActionSystem(reactions,rates)
	system.set_concentrations([2.,4.])
	system.run()
	print system.current_concentrations()
	print system.dydt()
	return

if __name__ == '__main__':
	main()