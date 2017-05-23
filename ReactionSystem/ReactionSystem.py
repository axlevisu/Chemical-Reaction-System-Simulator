# ReactionSystem.py
import numpy as np
from scipy.integrate import odeint


def odes(y,R,r,t=0):
	dydt = 
	return 


class ReactionSystem(object):
	"""docstring for ReactionSystem"""
	def __init__(self, reactions, rates):
		# super(ReactionSystem, self).__init__()
		self.reactions = reactions
		self.rates = rates

	def initialize(self,init_concentrations):
		if init_concentrations.shape[1] == self.reactions.shape[1]:
			self.init_concentrations = init_concentrations
			self.concentrations = init_concentrations
		else:
			print "Wrong initialization"

	def dydt(self):
		"""
		Returns rates of concentrations change at the current concentration
		"""
		return odes(self.concentrations,self.reactions,self.rates)

	def run(self,t=100000,ts=1000000):
		"""
		Run it for time t with ts time steps
		Outputs the concentration profile for the run
		default values are 100000 and 1000000
		"""	
		t = np.linspace(0, t, ts)
		y = self.concentrations
		output = 


