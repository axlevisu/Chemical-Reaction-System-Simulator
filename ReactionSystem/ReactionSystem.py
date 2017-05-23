# ReactionSystem.py
import numpy as np
from scipy.integrate import odeint

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def odes(y,t,R,k):
	dydt = np.zeros(y.shape[0])
	for i in xrange(R.shape[0]):
		l =  (R[i] > 0)*R[i]
		r = -((R[i] < 0)*R[i])
		# print l,r
		dydt += k[i]*(r-l)*arraypow(y,l)
	
	# print dydt
	return dydt 


class ReactionSystem(object):
	"""docstring for ReactionSystem"""
	def __init__(self, reactions, rates):
		# super(ReactionSystem, self).__init__()
		reactions = np.array(reactions)
		rates = np.array(rates)*1.0
		self.reactions = reactions
		self.rates = rates

#### TODO: Add a method to add more reactions
#### TODO: Add a method to visualize

	def set_concentrations(self,concentrations):
		concentrations = 1.0*np.array(concentrations)
		if concentrations.shape[0] == self.reactions.shape[1]:
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
		Run it for time t with ts time steps
		Outputs the concentration profile for the run
		default values are 100000 and 1000000
		"""	
		t = np.linspace(0, t, ts)
		y = self.concentrations
		output = odeint(odes, y, t, args= (self.reactions,self.rates))
		self.concentrations = output[-1,:]
		return output


def main():
	reactions = [[1,-1],[-1,1]]
	rates = [1,1]
	system = ReactionSystem(reactions,rates)
	system.set_concentrations([2,4])
	system.run()
	print system.current_concentrations()
	# print system.dydt()
	return

if __name__ == '__main__':
	main()